# Single-cell RNA-seq - normalization
# hbc-training 06_SC_SCT_normalization
# SJM 2024-07

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load previously cell+gene-filtered Seurat object, with joined and separate count layers in different seurat objects
load("data/cell_gene_filtered_seurat.RData")
load("data/cell_gene_filt_seurat_join.RData")

# Load cell cycle markers
load("data/cycle.rda")

# ---- Use of CellCycleScoring on a Seurat object in which the count matrices have been pre-filtered and reassigned causes mismatch errors -----#

# Before we make any comparisons across cells, we will apply a simple normalization. This is solely for the purpose of exploring the sources of variation in our data.
# Normalize the counts

seurat_norm <- NormalizeData(cell_gene_filt_seurat_join)
seurat_phase <- CellCycleScoring(seurat_norm,
                                g2m.features = g2m_genes,
                                s.features = s_genes
                                )
# seurat_norm <- CellCycleScoring(seurat_norm,
#                                  g2m.features = g2m_genes[which(g2m_genes %in% rownames(seurat_norm))],
#                                  s.features = s_genes[which(s_genes %in% rownames(seurat_norm))]
#                                 )

# The above leads to an error - presumably because there is a mismatch between the number of genes in the @features slot or elsewhere in the object compared to the underlying count matrices. which were filtered in the preceding QC script

# ----- #

# Load previously cell-filtered Seurat object (no gene filtering) to test if the gene filtering process in the QC script is causing issues with the @features slot
load("data/cell_filtered_seurat.RData")
cell_filt_seurat_join <- JoinLayers(cell_filtered_seurat)
seurat_norm <- NormalizeData(cell_filt_seurat_join)
seurat_phase <- CellCycleScoring(seurat_norm,
                                g2m.features = g2m_genes,
                                s.features = s_genes
)


### Note: importing two studies ('ctrl' and 'stim') as separate samples and layers within the same Seurat object may have been a reasonable strategy when these tutorials were written (July 2021) but appear to require splitting etc. with Seurat v5.0.0. Consider how best to keep each study as layers within a Seurat object as best practice in 2024, one corresponding to each experimental condition, with splitting as needed ater.

# ---- Splitting the Seurat object into separate objects for each sample set ----
# Split imported cell-gene-filtered Seurat object comprising two different sample sets back out into separate objects for easier downstream processing

split_seurat <- SplitObject(cell_gene_filtered_seurat, split.by = "sample")

# Score cells for cell cycle
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]])
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features = g2m_genes, s.features = s_genes)
}

# View cell cycle scores and phases assigned to cells
View(split_seurat[[1]]@meta.data)
View(split_seurat[[2]]@meta.data)

# In order to enable find variable selection, cell cycle phase commands below, do I now need to remerge the split objects back into a single object (seurat_phase)?

seurat_phase <- merge(split_seurat$ctrl, split_seurat$stim) # throws an error 2024-07-02
# -----

## Before performing PCA, this tutorial recommends selecting the most highly variable genes first and then scaling them so that the highest expressors don't unduly influence the dimension reduction

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")