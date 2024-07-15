# July/August 2021 workflow - Scot Matkovich June 2024
# HBC single-cell RNA-seq workshop
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/03_SC_quality_control-setup.md
# data at https://www.dropbox.com/s/vop78wq76h02a2f/single_cell_rnaseq.zip?dl=1

# Single-cell RNA-seq analysis - QC

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

# ---------------------
# Load data from a single sample manually using functions from Matrix and tidyverse packages

# Read in `matrix.mtx`
counts <- readMM("data/ctrl_raw_feature_bc_matrix/matrix.mtx.gz")

# Read in `genes.tsv`
genes <- read_tsv("data/ctrl_raw_feature_bc_matrix/features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cell_ids <- read_tsv("data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz", col_names = FALSE)$X1
# -----------------------

# Create a Seurat object for each sample
# The Read10X function reproduces the input procedures shown in the lines above, while the CreateSeuratObject function initiates a Seurat (S4?) object with slots
for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj) # change R object names on-the-fly
}

# Create a merged Seurat object
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, 
                       y = stim_raw_feature_bc_matrix, 
                       add.cell.id = c("ctrl", "stim"))

# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# calculate a 'novelty score'
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-") # note that this is a regular expression to search for markers in feature names. Change for non-human species or if mitochondrial genes are designated in a different way
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# -----
# Create metadata dataframe. Use this to make multiple changes to the metadata before incorporating back into the Seurat object
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

# Rename columns
metadata <- metadata |>
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
# -----

# Create .RData object to load at any time
# save(merged_seurat, file="data/merged_filtered_seurat.RData")

# -----
# Visualizations

# Visualize the number of cell counts per sample. From other information about this experiment, 12-13k cells are expected per sample
metadata |>
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# We see over 15,000 cells per sample, which is quite a bit more than the 12-13,000 expected. It is clear that we likely have some junk 'cells' present.

# Visualize the number UMIs/transcripts per cell
metadata |>
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
# The UMI counts per cell should generally be above 500, that is the low end of what we expect. If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.

# Visualize the distribution of genes detected per cell via histogram
metadata |> 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# We have similar expectations for gene detection as for UMI detection, although it may be a bit lower than UMIs. For high quality data, the proportional histogram should contain a single large peak that represents cells that were encapsulated. If we see a small shoulder to the left of the major peak (not present in our data), or a bimodal distribution of the cells, that can indicate a couple of things. It might be that there are a set of cells that failed for some reason. It could also be that there are biologically different types of cells (i.e. quiescent cell populations, less complex cells of interest), and/or one type is much smaller than the other (i.e. cells with high counts may be cells that are larger in size). Therefore, this threshold should be assessed with other metrics that we describe in this lesson.

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata |>
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
# We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score. The novelty score is computed by taking the ratio of nGenes over nUMI. If there are many captured transcripts (high nUMI) and a low number of genes detected in a cell, this likely means that you only captured a low number of genes and simply sequenced transcripts from those lower number of genes over and over again. These low complexity (low novelty) cells could represent a specific cell type (i.e. red blood cells which lack a typical transcriptome), or could be due to an artifact or contamination. Generally, we expect the novelty score to be above 0.80 for good quality cells.

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata |>
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark, unless of course you are expecting this in your sample.

# Visualize the correlation between genes detected and number of UMIs and determine whether there is a strong presence of cells with low numbers of genes/UMIs
metadata |>
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs (upper right quadrant of the plot). Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot. With this plot we also evaluate the slope of the line, and any scatter of data points in the bottom right hand quadrant of the plot. These cells have a high number of UMIs but only a few number of genes. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

# -----

# Filter out low quality cells using selected thresholds - these will change with experiment
cell_filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Create .RData object to load at any time
save(cell_filtered_seurat, file="data/cell_filtered_seurat.RData")

## Gene-level filtering
## Code below was originally for SeuratObject prior to 5.0.0, in which all counts were stored in one slot. The 5.0.0 and above version use 'layers' which need an alternate approach
# -----
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

# Reassign to filtered Seurat object
cell_gene_filt_seurat_join[['RNA']]$counts <- filtered_counts
cell_gene_filtered_seurat <- cell_gene_filt_seurat_join
cell_gene_filtered_seurat[['RNA']] <- split(cell_gene_filtered_seurat[['RNA']], f=cell_gene_filtered_seurat$sample)

# Create .RData objects to load at any time

save(cell_gene_filt_seurat_join, file="data/cell_gene_filt_seurat_join.RData")
save(cell_gene_filtered_seurat, file="data/cell_gene_filtered_seurat.RData")

# -----

# ## Gene-level filtering
# ## SeuratObject 5.0.0 and above
# # -----
# # Make copy of SeuratObject
# cell_gene_filtered_seurat <- cell_filtered_seurat
# # Get the names of all layers
# layer_names <- Layers(cell_gene_filtered_seurat)
# 
# # Loop over each layer in order to filter each layer separately
# for (layer in layer_names) {
#   # Extract the count matrix for this layer
#   counts <- LayerData(object = cell_gene_filtered_seurat, layer = layer)
#   
#   # Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
#   nonzero <- counts > 0
#   
#   # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
#   keep_genes <- Matrix::rowSums(nonzero) >= 10
#   
#   # Only keeping those genes expressed in more than 10 cells
#   filtered_counts <- counts[keep_genes, ]
#   
#   # Update the layer with the filtered count matrix
#   LayerData(object = cell_gene_filtered_seurat, layer = layer) <- filtered_counts
# }
# # -----
# 
# # -----
# # An alternative would be to join the layers, derive filtering parameters, and then apply the same set of filtering parameters to all layers
# cell_filtered_seurat_joined <- JoinLayers(cell_filtered_seurat)
# allcounts <- LayerData(object = cell_filtered_seurat_joined, layer = 'counts')
# rm(cell_filtered_seurat_joined)
# 
# # Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
# allnonzero <- allcounts > 0
# 
# # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
# keep_genes_all <- Matrix::rowSums(allnonzero) >= 10
# keep_genes_all_set <- names(which(keep_genes_all == T))
# 
# # Make copy of SeuratObject
# cell_gene_filtered_seurat <- cell_filtered_seurat
# # Get the names of all layers
# layer_names <- Layers(cell_gene_filtered_seurat)
# 
# # Loop over each layer in order to filter each layer separately, with the keep_genes_all gene set from above
# for (layer in layer_names) {
#   # Extract the count matrix for this layer
#   counts <- LayerData(object = cell_gene_filtered_seurat, layer = layer)
#   
#   # Only keeping those genes expressed in more than 10 cells from prior filtering workflow
#   filtered_counts <- counts[which(rownames(counts) %in% keep_genes_all_set), ]
#   
#   # Update the layer with the filtered count matrix
#   LayerData(object = cell_gene_filtered_seurat, layer = layer) <- filtered_counts
# }
# # -----

## Exercises
# Save (cell)-filtered subset to new metadata
metadata_clean <- cell_gene_filtered_seurat@meta.data

# Perform previous QC ggplot steps using metadata_clean in place of metadata


