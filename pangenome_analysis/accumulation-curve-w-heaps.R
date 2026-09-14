library("micropan")
library(dplyr)
library(ggplot2)
set.seed(123)
global_filename <- "shared_pangenome_matrix.csv"
british_filename <- "british_isles.csv"
# Read files
british_meta <- read.csv(british_filename)
pa <- read.csv(global_filename, header=T, row.names = 1)
# Subset BI pangenome
british_total <- pa[,colnames(pa) %in% british_meta$new.ID]
british_total <- british_total[rowSums(british_total==0, na.rm =TRUE)<ncol(british_total),]

# Store final data
bi_pa <- t(british_total)
bi_pa_nonsing <- bi_pa[,colSums(bi_pa) != 1]
pa <- t(pa)
pa_nonsing <- pa[,colSums(pa) != 1]


get_accumulation_data <- function(pg_matrix, identifier, permutations=100){
  # Generate rarefaction data
  line_data <- rarefaction(pg_matrix, n.perm = permutations)
  heap_data <- heaps(pg_matrix, n.perm = permutations)
  heap=unname(signif(heap_data[2],2))
  heap_out <- cbind(identifier,heap)
  # Identify the median values across the permutations
  perm_data <- line_data[2:nrow(line_data),2:ncol(line_data)]
  # Get median data
  mean_perm <- round(apply(perm_data, 1, mean, na.rm=T))
  # Get minimum data
  min_perm <- round(apply(perm_data, 1, min, na.rm=T))
  # Get maximum data
  max_perm <- round(apply(perm_data, 1, max, na.rm=T))
  # Combine
  mean_data <- as.data.frame(cbind(1:length(mean_perm), mean_perm, min_perm, max_perm, rep(identifier, length(mean_perm))))
  return(list("mean_data"=mean_data, "heap_data"=heap_out))
}

global_acc <- get_accumulation_data(pa, "Global", permutations = 1000)
global_nonsing_acc <- get_accumulation_data(pa_nonsing, "Global (w/o Singletons)", permutations = 1000)
bi_acc <- get_accumulation_data(bi_pa, "British Isles", permutations = 1000)
bi_nonsing_acc <- get_accumulation_data(bi_pa_nonsing, "British Isles (w/o Singletons)", permutations = 1000)

total_heap_data <- rbind(global_acc$heap_data, global_nonsing_acc$heap_data, bi_acc$heap_data, bi_nonsing_acc$heap_data)



total_data <- rbind(global_acc$mean_data, global_nonsing_acc$mean_data, bi_acc$mean_data, bi_nonsing_acc$mean_data)
total_data[,1:4] <- mutate_all(total_data[,1:4], function(x) as.numeric(as.character(x)))
total_data$V5 <- factor(total_data$V5, levels=c("Global", "British Isles", "Global (w/o Singletons)", "British Isles (w/o Singletons)"))
p <- ggplot(total_data) +
  geom_ribbon(aes(x=V1, ymin=min_perm, ymax=max_perm, fill = V5), alpha=0.5) +
  geom_line(aes(x=V1, y=mean_perm, color = V5)) +
  xlim(c(0,1200)) +
  xlab("Number of genomes") +
  ylab("Number of gene families") +
  scale_color_manual(name = "Dataset",
                      labels = c("Global", "British Isles", "Global (w/o Singletons)", "British Isles (w/o Singletons)"),
                      values = c("#E31733", "#031C84", "#EC2190", "#006682")) +   
  scale_fill_manual(name = "Dataset",
                    labels = c("Global", "British Isles", "Global (w/o Singletons)", "British Isles (w/o Singletons)"),
                    values = c("#FBB1BB", "#319AFE", "#E487C8", "#00D9DB")) +  
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        # explicitly set the horizontal lines (or they will disappear too)
        legend.title = element_text(face="bold", color = "black", size = 14),
        legend.text = element_text(face="bold", color = "black", size = 12),
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=14, colour = "black"),
        axis.title = element_text(face = "bold", size = 16, colour = "black"))

output_filename <- "/Users/user/Documents/manchester/global/results/gene-family-pangenome/pangenome_analysis/gene_family.acc_curve.pdf"
pdf(file =output_filename,width = 10, height = 6)
print(p)
dev.off()

total_heap_data <- as.data.frame(total_heap_data)
total_heap_data[,2] <- as.numeric(as.character(total_heap_data[,2]))
total_heap_data$identifier <- factor(total_heap_data$identifier, levels=c("Global", "British Isles", "Global (w/o Singletons)", "British Isles (w/o Singletons)"))
bp <- ggplot(data=total_heap_data) +
  geom_bar(stat="identity", aes(x=identifier, y=heap, fill=identifier)) +
  ylim(c(0,max(total_heap_data$heap) + 0.1)) +
  xlab("Dataset") +
  ylab("Heap value") +
  scale_fill_manual(name = "Dataset",
                     labels = c("Global", "British Isles", "Global (w/o Singletons)", "British Isles (w/o Singletons)"),
                     values = c("#E31733", "#031C84", "#EC2190", "#006682")) + 
  scale_x_discrete(labels=c("G", "BI", "G (w/o S)", "BI (w/o S)")) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        # explicitly set the horizontal lines (or they will disappear too)
        legend.title = element_text(face="bold", color = "black", size = 14),
        legend.text = element_text(face="bold", color = "black", size = 12),
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=14, colour = "black"),
        axis.title = element_text(face = "bold", size = 16, colour = "black"))

output_filename <- "/Users/user/Documents/manchester/global/results/gene-family-pangenome/pangenome_analysis/gene_family.heap.pdf"
pdf(file =output_filename,width = 10, height = 6)
print(bp)
dev.off()

save.image(file = "/Users/user/Documents/manchester/global/results/gene-family-pangenome/pangenome_analysis/acc_results.RData")


