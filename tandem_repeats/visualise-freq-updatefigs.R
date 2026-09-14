library(ggplot2)
library(data.table)

# Load colorscheme
colourscheme_filename <- "/Users/user/Documents/manchester/global/data/metadata/colourscheme.csv"
colourscheme <- read.csv(colourscheme_filename, header=T)
cluster_colourscheme <- colourscheme[colourscheme$type == 'cluster_id',]
rownames(cluster_colourscheme) <- as.character(cluster_colourscheme$value)
cluster_col_dict <- as.list(cluster_colourscheme$colour)
# Character colorscheme
char_col_dict = list("Cluster 1" = "#340043",
                     "Cluster 2" = "#2D3E78",
                     "Cluster 3" =  "#1E7F79",
                     "Cluster 4" = "#4FC150",
                     "Cluster 5" = "#FCE51D")


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


metadata <- read.csv("~/Documents/manchester/global/data/metadata/microreact.shared.csv", header=F, row.names=1)

repeat_unit_data <- read.csv("/Users/user/Documents/manchester/tr/results/overview/consensus_size_freq.csv")

total_tr_data <- read.csv("/Users/user/Documents/manchester/tr/results/overview/total_tandem_repeat.csv")

repeat_length_data <- read.csv("/Users/user/Documents/manchester/tr/results/overview/repeat_length_freq.csv")


total_tr_data$Repeat.Label <- rep("Short", nrow(total_tr_data))
total_tr_data$Repeat.Label[total_tr_data$Consensus.Size >= 10] <- "Medium"
total_tr_data$Repeat.Label[total_tr_data$Consensus.Size >= 50] <- "Long"

total_tr_data$Length.Label <- rep("Short", nrow(total_tr_data))
total_tr_data$Length.Label[total_tr_data$Consensus.Size >= 30] <- "Medium"
total_tr_data$Length.Label[total_tr_data$Consensus.Size >= 100] <- "Long"

# No. TRs overall
total_tr_count <- as.data.frame(table(total_tr_data$Isolate))
# Add cluster info
cluster <- c()
for(iso in total_tr_count$Var1){
  cluster <- c(cluster,
               unique(total_tr_data$Cluster[total_tr_data$Isolate == iso]))
}
total_tr_count$Cluster <- factor(cluster)
# Add genotype info
genotype <- c()
for(iso in total_tr_count$Var1){
  genotype <- c(genotype, metadata[iso,]$V7)
}
total_tr_count$Genotype <- factor(genotype)

total_tr_P <- ggplot(total_tr_count, aes(x=Cluster, y=Freq, fill=Cluster)) + 
  geom_violin(outlier.shape = NA, linewidth = 1, alpha=0.5) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(limits = quantile(total_tr_count$Freq, c(0.1, 0.9))) +
  labs(x = "Cluster", y = "Total No. TR's") +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  #scale_colour_manual(values=unlist(char_col_dict)) + 
  scale_fill_manual(values=char_col_dict) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))
output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/total_TR_count.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(total_tr_P)
dev.off()


total_tr_G_P <- ggplot(total_tr_count, aes(x=Freq, y=Genotype, fill=Genotype)) + 
  geom_violin(outlier.shape = NA, linewidth = 1, alpha=0.5) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  scale_x_continuous(limits = quantile(total_tr_count$Freq, c(0.1, 0.9))) +
  labs(x = "Total No. TR's", y = "Genotype") +
  #scale_colour_manual(values=unlist(char_col_dict)) + 
  scale_fill_manual(values=c25) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))
output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/total_TR_count-genotype.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(total_tr_G_P)
dev.off()

label_tr_data <- total_tr_data[, c("Isolate", "Repeat.Label")]
label_frequency <- as.data.frame(table(label_tr_data))

length_label_tr_data <- total_tr_data[, c("Isolate", "Length.Label")]
length_label_frequency <- as.data.frame(table(length_label_tr_data))
# Add cluster information
cluster <- c()
for(iso in label_frequency$Isolate){
  cluster <- c(cluster,
               unique(repeat_unit_data$Cluster[repeat_unit_data$Isolate == iso]))
}
label_frequency$Cluster <- cluster


cluster <- c()
for(iso in length_label_frequency$Isolate){
  cluster <- c(cluster,
               unique(repeat_unit_data$Cluster[repeat_unit_data$Isolate == iso]))
}
length_label_frequency$Cluster <- cluster



average_unit_size = list("Isolate" = c(),
                         "Mean" = c(),
                         "Cluster" = c())
average_repeat_length = list("Isolate" = c(),
                         "Mean" = c(),
                         "Cluster" = c())
average_repeat_class = list("Isolate" = c(),
                            "Label" = c(),
                             "Mean" = c(),
                             "Cluster" = c())
average_num_repeats = list("Isolate" = c(),
                             "Mean" = c(),
                             "Cluster" = c())
total_size_repeats = list("Isolate" = c(),
                           "Total" = c(),
                           "Cluster" = c())
for(iso in total_tr_count$Var1){
  avg_unit_size <- mean(total_tr_data$Consensus.Size[total_tr_data$Isolate == iso])
  avg_rep_len <- mean(total_tr_data$Repeat.Length[total_tr_data$Isolate == iso])
  avg_n_rep <- mean(total_tr_data$Copies[total_tr_data$Isolate == iso])
  
  total_tr_len <- sum(total_tr_data$Repeat.Length[total_tr_data$Isolate == iso])
  avg_unit_size <- mean(total_tr_data$Consensus.Size[total_tr_data$Isolate == iso])
  
  
  
  cluster <-  unique(total_tr_data$Cluster[total_tr_data$Isolate == iso])
  
  average_unit_size$Isolate <- c(average_unit_size$Isolate, iso)
  average_unit_size$Mean <- c(average_unit_size$Mean, avg_unit_size)
  average_unit_size$Cluster <- c(average_unit_size$Cluster, cluster)
  
  average_repeat_length$Isolate <- c(average_repeat_length$Isolate, iso)
  average_repeat_length$Mean <- c(average_repeat_length$Mean, avg_rep_len)
  average_repeat_length$Cluster <- c(average_repeat_length$Cluster, cluster)
  
  average_num_repeats$Isolate <- c(average_num_repeats$Isolate, iso)
  average_num_repeats$Mean <- c(average_num_repeats$Mean, avg_n_rep)
  average_num_repeats$Cluster <- c(average_num_repeats$Cluster, cluster)
  
  total_size_repeats$Isolate <- c(total_size_repeats$Isolate, iso)
  total_size_repeats$Total <- c(total_size_repeats$Total, total_tr_len)
  total_size_repeats$Cluster <- c(total_size_repeats$Cluster, cluster)  
  }

average_num_rep_df <- as.data.frame(average_num_repeats)
average_rep_length_df <- as.data.frame(average_repeat_length)
average_unit_size_df <- as.data.frame(average_unit_size)
total_size_df <- as.data.frame(total_size_repeats)

avg_rep_length_P <- ggplot(average_rep_length_df, aes(x=Cluster, y=Mean, fill=Cluster)) + 
  geom_violin(outlier.shape = NA, linewidth = 1, alpha=0.5) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(limits = quantile(average_rep_length_df$Mean, c(0.1, 0.9))) +
  labs(x = "Cluster", y = "Average Repeat Length") +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  #scale_colour_manual(values=unlist(char_col_dict)) + 
  scale_fill_manual(values=char_col_dict) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))

output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/average_repeat_length.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(avg_rep_length_P)
dev.off()



average_num_rep_P <- ggplot(average_num_rep_df, aes(x=Cluster, y=Mean, fill=Cluster)) + 
  geom_violin(outlier.shape = NA, linewidth = 1, alpha=0.5) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(limits = quantile(average_num_rep_df$Mean, c(0.1, 0.9))) +
  labs(x = "Cluster", y = "Average Number of Repeat Copies") +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  #scale_colour_manual(values=unlist(char_col_dict)) + 
  scale_fill_manual(values=char_col_dict) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))
output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/average_number_of_repeats.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(average_num_rep_P)
dev.off()



average_unit_size_P <- ggplot(average_unit_size_df, aes(x=Cluster, y=Mean, fill=Cluster)) + 
  geom_violin(outlier.shape = NA, linewidth = 1, alpha=0.5) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  scale_y_continuous(limits = quantile(average_unit_size_df$Mean, c(0.1, 0.9))) +
  labs(x = "Cluster", y = "Average Repeat Unit Size") +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  #scale_colour_manual(values=unlist(char_col_dict)) + 
  scale_fill_manual(values=char_col_dict) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))

output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/average_unit_size.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(average_unit_size_P)
dev.off()

size_label_count_P <- ggplot(label_frequency, aes(y=Freq, fill=Repeat.Label)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(label_frequency$Freq, c(0.1, 0.9))) +
  facet_grid(. ~ Cluster) +
  ylab("Count") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))
  
output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/count_by_size_category.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(size_label_count_P)
dev.off()


length_label_count_P <- ggplot(length_label_frequency, aes(y=Freq, fill=Length.Label)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(length_label_frequency$Freq, c(0.1, 0.9))) +
  facet_grid(. ~ Cluster) +
  ylab("Count") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        legend.position = "none",
        panel.grid.major.y = element_line( size=0.1, color="black", linetype="dashed" ) ,
        axis.line = element_line(colour = "black"),
        axis.text = element_text(face="bold", size=5, colour = "black"),
        axis.title = element_text(face = "bold", size = 6, colour = "black"))
output_filename <- "/Users/user/Documents/manchester/tr/figures/figure_update/count_by_length_category.pdf"
pdf(file=output_filename,width = 3.14, height = 3.14)
print(length_label_count_P)
dev.off()

