# Single-cell RNA-seq - normalization
# hbc-training 06_SC_SCT_normalization
# SJM 2024-07

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load cell cycle markers
load("data/cycle.rda")

# Load previously cell-filtered Seurat object (no gene filtering) to test if the gene filtering process in the QC script is causing issues with the @features slot
load("data/cell_filtered_seurat.RData")
cell_filt_seurat_join <- JoinLayers(cell_filtered_seurat)
seurat_norm <- NormalizeData(cell_filt_seurat_join)
seurat_phase <- CellCycleScoring(seurat_norm,
                                g2m.features = g2m_genes,
                                s.features = s_genes
)

### Note: importing two studies ('ctrl' and 'stim') as separate samples and layers within the same Seurat object may have been a reasonable strategy when these tutorials were written (July 2021) but appear to require splitting etc. with Seurat v5.0.0. Consider how best to keep each study as layers within a Seurat object as best practice in 2024, one corresponding to each experimental condition, with splitting as needed later.

# ----
# Alternative for CellCycleScoring followed by gene filtering
# 
# Perform CellCycleScoring as above on the cell-filtered object. Then, once these designations are present in the @meta.data slot, perform gene-filtering as in the QC script
#
# Gene-level filtering
# Code below was originally for SeuratObject prior to 5.0.0, in which all counts were stored in one slot. The 5.0.0 and above version use 'layers' which need an alternate approach
# 
# Join layers
cell_gene_filt_seurat_join <- JoinLayers(cell_filtered_seurat)
# Extract counts
counts <- GetAssayData(object = cell_gene_filt_seurat_join, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to desired Seurat object
seurat_phase[['RNA']]$counts <- filtered_counts
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

# What about mitochondrial fraction rather than cell cycle phase?

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)
# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))
# Plot the PCA colored by mitoFraction
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")
################################################################################

# ---- Use of CellCycleScoring on a Seurat object in which the count matrices have been pre-filtered and reassigned causes mismatch errors -----#

# Load previously cell+gene-filtered Seurat object, with joined and separate count layers in different seurat objects
# load("data/cell_gene_filtered_seurat.RData")
# load("data/cell_gene_filt_seurat_join.RData")

# Before we make any comparisons across cells, we will apply a simple normalization. This is solely for the purpose of exploring the sources of variation in our data.
# Normalize the counts

# seurat_norm <- NormalizeData(cell_gene_filt_seurat_join)
# seurat_phase <- CellCycleScoring(seurat_norm,
#                                 g2m.features = g2m_genes,
#                                 s.features = s_genes
#                                 )

# seurat_norm <- CellCycleScoring(seurat_norm,
#                                  g2m.features = g2m_genes[which(g2m_genes %in% rownames(seurat_norm))],
#                                  s.features = s_genes[which(s_genes %in% rownames(seurat_norm))]
#                                 )

# The above leads to an error - presumably because there is a mismatch between the number of genes in the @features slot or elsewhere in the object compared to the underlying count matrices. which were filtered in the preceding QC script

# ----- #
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

# Or in a workflow like this, should highly variable genes and PCA be performed for each list element (Seurat object) separately, if they can't be remerged?
# -----

################################################################################

# Normalization and regressing out sources of unwanted variation using SCTransform
# Now that we have established which effects are observed in our data, we can use the SCTransform method to regress out these effects. The SCTransform method was proposed as a better alternative to the log transform normalization method that we used for exploring sources of unwanted variation. The method not only normalizes data, but it also performs a variance stabilization and allows for additional covariates to be regressed out.
# 
# As described earlier, all genes cannot be treated the same. As such, the SCTransform method constructs a generalized linear model (GLM) for each gene with UMI counts as the response and sequencing depth as the explanatory variable. Information is pooled across genes with similar abundances, to regularize parameter estimates and obtain residuals which represent effectively normalized data values which are no longer correlated with sequencing depth

# Iterating over samples in a dataset
# Since we have two samples in our dataset (from two conditions), we want to keep them as separate objects and transform them as that is what is required for integration. We will first split the cells in seurat_phase object into "Control" and "Stimulated":

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# change allowable R object size from 500 Mb to 4000 Mb
options(future.globals.maxSize = 4000 * 1024^2)

# the SCTransform function below carries out the GLM above (UMI counts as response, sequencing depth as explanatory variable) with the addition of mitoRatio as a further variable to regress out
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
}

# View summmaries
split_seurat$ctrl@assays
split_seurat$stim@assays

# Save the split seurat object (load back later with loadRDS)
saveRDS(split_seurat, "data/split_seurat_phase_mitoRatio.rds")
