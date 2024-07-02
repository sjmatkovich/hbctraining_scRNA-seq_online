# Single-cell RNA-seq - normalization
# hbc-training 06_SC_SCT_normalization
# SJM 2024-07

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load previously cell+gene-filtered Seurat object
load("data/seurat_cell_gene_filtered.RData")

# Before we make any comparisons across cells, we will apply a simple normalization. This is solely for the purpose of exploring the sources of variation in our data.

# Normalize the counts
seurat_phase <- NormalizeData(cell_gene_filtered_seurat)

# Load cell cycle markers
load("data/cycle.rda")

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# Score cells for cell cycle
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features = g2m_genes, s.features = s_genes)
}

# View cell cycle scores and phases assigned to cells                             
View(split_seurat[[1]]@meta.data)                             
View(split_seurat[[2]]@meta.data)

