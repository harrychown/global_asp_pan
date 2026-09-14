# Analyse pan-genomic data
library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(plyr)
library(tidyverse)
library(reshape2)
set.seed(1234)
# Set input
exonerate_filename <- "~/Downloads/panstripe/reference-matrix.csv"
orthogroup_filename <- "~/Downloads/panstripe/Orthogroups.tsv"
singleton_filename <- "~/Downloads/panstripe/Orthogroups_UnassignedGenes.tsv"
phylogeny_filename <- "~/Downloads/panstripe/RAxML_result.globalpaper_MP"
metadata_filename <- "~/Downloads/panstripe//globalMicroreact.csv"
colourscheme_filename <- "~/Downloads/panstripe/colourscheme.csv"


# Read input
exonerate <- read.csv(exonerate_filename, header=T, row.names = 1)
orthogroup <- read.table(orthogroup_filename, sep ='\t', header=T, row.names = 1, na.strings="")
singleton <- read.table(singleton_filename, sep ='\t', header=T, row.names = 1, na.strings="")
colourscheme <- read.csv(colourscheme_filename, header=T)
metadata <- read.csv(metadata_filename, header=T)
tree <- read.tree(phylogeny_filename)
# Re-name and re-order columns
rename_col <- function(x){
  col <- colnames(x)
  col <- sapply(strsplit(col, "\\."), '[', 1)
  colnames(x) <-  col
  x <- x[,sort(col)]
  return(x)
}

# Binarise data
binarise <- function(x){
  x[is.na(x)] <- as.numeric(0)
  x[x != 0] <- as.numeric(1)
  x[] <- lapply(x, as.numeric)
  return(x)
}
# Add 'dummy' AF293 column to exonerate
exonerate$AF293 <- rep(1, nrow(exonerate))
# Re-order
exonerate <- rename_col(exonerate)
orthogroup <- rename_col(orthogroup)
singleton <- rename_col(singleton)

# Combine matrices
bin_exonerate <- binarise(exonerate)
bin_orthogroup <- binarise(orthogroup)
bin_singleton <- binarise(singleton)
total <- rbind(bin_exonerate, bin_orthogroup, bin_singleton)

# Match tree with matrix
shared_tree_data <- tree$tip.label[tree$tip.label %in% colnames(total)]
shared_total <- total[,shared_tree_data]
shared_tree <-drop.tip(tree,tree$tip.label[-match(shared_tree_data, tree$tip.label)])
write.tree(shared_tree, "shared_tree.nwk")
write.csv(shared_total, "shared_pangenome.csv", quote = F)

