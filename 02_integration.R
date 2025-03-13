# 02_integration.R
#
# This script performs normalization and integration of the scRNA-seq data:
# 1. SCTransform normalization for each sample
# 2. Selection of integration features
# 3. Harmony integration to correct for batch effects
# 4. PCA and UMAP dimensionality reduction

# Load required libraries
library(Seurat)
library(harmony)
library(future)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Enable parallel processing
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)

# Load preprocessed data
cat("Loading preprocessed data...\n")
pd <- readRDS("data/merged_seurat_object_preprocessed.rds")

# Split by sample for integration
cat("Splitting object by sample ID for integration...\n")
pd_split <- SplitObject(pd, split.by = "SamplID")

# SCTransform normalization on each sample
cat("Performing SCTransform normalization on each sample...\n")
for (i in names(pd_split)) {
  cat(paste0("Processing sample: ", i, "\n"))
  pd_split[[i]] <- SCTransform(pd_split[[i]], verbose = TRUE, vars.to.regress = "percent.mt")
}

# Select integration features
cat("Selecting integration features...\n")
features <- SelectIntegrationFeatures(object.list = pd_split)

# Prepare for integration
cat("Preparing for SCT integration...\n")
pd_split <- PrepSCTIntegration(object.list = pd_split, anchor.features = features, verbose = FALSE)

# Merge the split objects back together
cat("Merging objects for Harmony integration...\n")
pd_merged <- merge(pd_split[[1]], 
                   y = pd_split[2:length(pd_split)], 
                   project = "PD_scRNA", 
                   merge.data = TRUE)

# Run PCA for Harmony integration
cat("Running PCA...\n")
pd_merged <- RunPCA(object = pd_merged, npcs = 20, assay = "SCT", features = features)

# Visualize PCA results
cat("Visualizing PCA results...\n")
pdf("figures/pca_heatmap.pdf", width = 12, height = 10)
DimHeatmap(pd_merged, dims = 1:15, cells = 1000, balanced = TRUE)
dev.off()

# Run Harmony integration
cat("Running Harmony integration...\n")
pd_integrated <- RunHarmony(object = pd_merged,
                            assay.use = "SCT",
                            reduction = "pca",
                            dims.use = 1:15,
                            group.by.vars = "SamplID",
                            plot_convergence = TRUE)

# Run UMAP on Harmony-integrated data
cat("Running UMAP on Harmony-integrated data...\n")
pd_integrated <- RunUMAP(object = pd_integrated, assay = "SCT", reduction = "harmony", dims = 1:15)

# Build neighbor graph for clustering
cat("Building neighbor graph...\n")
pd_integrated <- FindNeighbors(object = pd_integrated, assay = "SCT", reduction = "harmony", graph.name = "harmony_snn", dims = 1:15)

# Visualize integration
cat("Creating visualization plots...\n")
p1 <- DimPlot(pd_integrated, group.by = "Mutation", reduction = "umap") + ggtitle("Integration by Mutation")
p2 <- DimPlot(pd_integrated, group.by = "Region", reduction = "umap") + ggtitle("Integration by Region")
p3 <- DimPlot(pd_integrated, group.by = "SamplID", reduction = "umap") + ggtitle("Integration by Sample")

# Check feature distribution
p4 <- FeaturePlot(pd_integrated, reduction = 'umap', features = c("nFeature_RNA", "percent.mt"))

# Save plots
ggsave("figures/integration_by_mutation.pdf", plot = p1, width = 10, height = 8)
ggsave("figures/integration_by_region.pdf", plot = p2, width = 10, height = 8)
ggsave("figures/integration_by_sample.pdf", plot = p3, width = 10, height = 8)
ggsave("figures/qc_post_integration.pdf", plot = p4, width = 10, height = 8)

# Save the integrated object
cat("Saving integrated Seurat object...\n")
saveRDS(pd_integrated, file = "data/integrated_seurat_object.rds")

cat("Integration completed! The integrated Seurat object is saved at data/integrated_seurat_object.rds\n")