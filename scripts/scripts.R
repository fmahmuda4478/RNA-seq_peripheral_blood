# Install GEOquery from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("GEOquery")

# Load the package
library(GEOquery)

# Example: Download a GEO Series dataset
gse <- getGEO("GSE21942", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Check the structure
str(gse)
# Extract expression data
expr_data <- exprs(gse[[1]])  # Use [[1]] for the first platform

# View sample information
sample_info <- pData(gse[[1]])

# View feature annotation
feature_info <- fData(gse[[1]])

# Summary of the dataset
summary(expr_data)

# Download raw files
getGEOSuppFiles("GSE21942")


save(gse, file = "GSE2)



gse[["GSE21942_series_matrix.txt.gz"]]

geo_data<- getGEO("GSE21942",GSEMatrix = TRUE,getGPL =TRUE)

# Basic heatmap of the expression data
heatmap(expr_data[1:100, ], main = "Heatmap of Top 100 Genes")

save(gse, file = "GSE21942.RData")

# Load it back
load("GSE21942.RData")

install.packages("tidyverse")

install.packages("ggpubr")

install.packages("openxlsx")
 install.packages("BiocManager")


BiocManager::install("TCGAiolinks")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# Install DESeq2
BiocManager::install("DESeq2")

# Load the package
library(DESeq2)

# Example data
library(DESeq2)
countData <- matrix(rpois(1000, lambda = 10), ncol = 5)
colData <- data.frame(condition = factor(c("A", "A", "B", "B", "A")), row.names = colnames(countData))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

# Estimate size factors (normalization)
dds <- estimateSizeFactors(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Variance stabilizing transformation
vsd <- vst(dds)

# or use rlog(dds) for regularized log transformation
normalized_matrix <- assay(vsd)

# Example raw count matrix
counts <- matrix(rpois(1000, lambda = 5), nrow = 200, ncol = 5)  # 200 genes x 5 samples
rownames(counts) <- paste0("Gene", 1:200)

# Filter: Keep genes with a minimum count threshold in at least N samples
min_count <- 10  # Minimum count threshold
min_samples <- 2  # Number of samples that must meet the threshold

filtered_counts <- counts[rowSums(counts >= min_count) >= min_samples, ]

# View dimensions before and after filtering
dim(counts)  # Original dimensions
dim(filtered_counts)  # Dimensions after filtering

library(DESeq2)

# Example data
countData <- matrix(rpois(1000, lambda = 5), nrow = 200, ncol = 5)
rownames(countData) <- paste0("Gene", 1:200)
colData <- data.frame(condition = factor(c("A", "A", "B", "B", "A")))

# Filter: Keep genes with a row mean above a threshold
filtered_countData <- countData[rowMeans(countData) > 1, ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_countData, colData = colData, design = ~ condition)


install.packages("ggplot2")
library(ggplot2)
countsData |>
   rename(gene=x)

ggplot(data=,aes(x=,y=)) +geom_type()
plot_data <- countData %>%
  select(Gene_Name = Gene_Symbol, Sample_1 = GSM12345, Sample_2 = GSM12346) %>% 
  pivot_longer(cols = c(Sample_1, Sample_2), names_to = "Sample", values_to = "Count")





