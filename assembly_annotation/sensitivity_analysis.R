library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(retroPal)
library(nortest)

b90_k90_filename <- "/Users/user/Documents/manchester/global/submission/rebuttal/sensitivity_analysis/b90_k90/b90_k90.txt"
b95_k95_filename <- "/Users/user/Documents/manchester/global/submission/rebuttal/sensitivity_analysis/b95_k95/b95_k95.txt"
pangenome_file <- "/Users/user/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Documents/manchester/pangenome_evolution/data/pangenome/orthofinder/shared_orthofinder_matrix.csv"
asm_stats_file <- "~/Documents/manchester/global/submission/rebuttal/filtered.assembly-stats.csv"
output_folder <- "/Users/user/Documents/manchester/global/submission/rebuttal/sensitivity_analysis"
# Read in files
metadata <- read.csv("/Users/user/Documents/manchester/global/data/metadata/globalMicroreact.csv")
pg <- read.csv(pangenome_file, row.names = 1)
b90_k90 <- read.table(b90_k90_filename)
b95_k95 <- read.table(b95_k95_filename)
asm_stats <- read.csv(asm_stats_file)



# Example column names — adjust to yours:
# scaffold_len, n_genes, n_denovo_removed
asm_stats$denovo_pre_filter <- asm_stats$denovo_genes_dropped + asm_stats$denovo_gene_count
asm_stats$denovo_pre_filter_log <- log10(asm_stats$denovo_pre_filter + 1)
asm_stats$denovo_gene_count_log <- log10(asm_stats$denovo_gene_count + 1)

df_long <- asm_stats %>%
  mutate(scaffold_len = as.numeric(scaffold_sequence_total)) %>%
  pivot_longer(
    cols = c(denovo_pre_filter, denovo_gene_count),
    names_to  = "metric",
    values_to = "count"
  ) %>%
  mutate(
    metric = recode(metric,
                    denovo_pre_filter = "De novo genes, pre-filter",
                    denovo_gene_count = "De novo genes, post-filter")
  )

df_long$metric <- factor(df_long$metric, levels = c("De novo genes, pre-filter", "De novo genes, post-filter"))

stats <- df_long %>%
  group_by(metric) %>%
  do({
    m <- lm(count ~ scaffold_len, data = .)
    s <- summary(m)
    tibble(
      r2 = s$r.squared,
      p  = coef(s)[2, "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(label = sprintf("R² = %.3f\np = %.2g", r2, p))


gene_hist <- ggplot(df_long, aes(x = count)) +
  geom_histogram(bins = 50, colour = "black", fill = "grey70") +
  facet_wrap(~ metric, scales = "free_x") +
  labs(x = "Count", y = "Frequency") +
  pub_theme(y_axis_labels=T, x_axis_labels = T)

gene_QQ <- ggplot(df_long, aes(sample = count)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ metric, scales = "free") +
  xlab("Theoretical Quantiles") +
  ylab("Sampel Quantiles") +
  pub_theme(y_axis_labels=T, x_axis_labels = T)


gene_count <- ggplot(df_long, aes(x = scaffold_sequence_total, y=count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ metric, scales = "free_x") +
  geom_text(
    data = stats,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = 5*.3527
  ) +
  xlab("Assembly Length (bp)") +
  ylab("Count") +
  pub_theme(y_axis_labels=T, x_axis_labels = T)

ggsave(file=file.path(output_folder, "gene_QQ.svg"), plot=gene_QQ, width=8, height=8, units="cm", bg="white")
ggsave(file=file.path(output_folder, "gene_count.svg"), plot=gene_count, width=8, height=8, units="cm", bg="white")



# Subset based on stringency
b80_k80_meta <- metadata[which(metadata$ID %in% colnames(pg)),]
b90_k90_meta <- metadata[which(metadata$ID %in% b90_k90[,1]),]
b95_k95_meta <- metadata[which(metadata$ID %in% b95_k95[,1]),]

# Display geographical distribution between samples
continent_distributions <- list( b80_k80 = table(b80_k80_meta$Continent),
      b90_k90 = table(b90_k90_meta$Continent),
      b95_k95 = table(b95_k95_meta$Continent))

cluster_distributions <- list( b80_k80 = table(b80_k80_meta$Cluster),
                                 b90_k90 = table(b90_k90_meta$Cluster),
                                 b95_k95 = table(b95_k95_meta$Cluster))

as.data.frame(cluster_distributions)
as.data.frame(continent_distributions)[,c(2,4,6)]

# Calculate pangenome metrics with sub-samples
b80_k80_pg <- pg
b90_k90_pg <- pg[,which(b90_k90[,1] %in% colnames(pg))]
b95_k95_pg <- pg[,which(b95_k95[,1] %in% colnames(pg))]

b90_k90_pg <- b90_k90_pg[rowSums(b90_k90_pg) != 0,]
b95_k95_pg <- b95_k95_pg[rowSums(b95_k95_pg) != 0,]

