# Script to analyse the coverages of genes (per strain)
args = commandArgs(trailingOnly=TRUE)
input_filename <- args[1]
output_filename <-  args[2]
upper_outlier_filename <- args[3]
raw <- read.table(input_filename)
# Extract reference mitochondrial encoded genes
mitogene <- raw[grep("AfuMt", raw[,1]),]
# Extract reference nuclear encoded genes
afu_generic <- raw[grep("Afu", raw[,1]),1]
nucl_id <- afu_generic[grep("AfuMt", afu_generic, invert=TRUE)]
nuclgene <- raw[which(raw[,1] %in% nucl_id),]
# Extract de novo annotated genes
brakergene <- raw[grep("Afu", raw[,1], invert=T),]
# Calculate threshold for outliers
Q1 <- as.numeric(unname(summary(nuclgene[,2])[2]))
IQR <-  IQR(nuclgene[,2])
limit <- Q1 - (1.5 * IQR)
outliers <- brakergene[which(brakergene[,2] < limit),1]
if(length(outliers>0)){
  write(outliers, output_filename)
}
# Identify which mitochondrial genes are upper bound outliers
Q3 <- as.numeric(unname(summary(nuclgene[,2])[5]))
upper_limit <- Q3 + (1.5 * IQR)
upper_MT <- mitogene[which(mitogene[,2] > limit),1]
if(length(upper_MT>0)){
  write(upper_MT, upper_outlier_filename)
}

