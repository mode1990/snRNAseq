# 03_clustering.R
#
# This script performs clustering and cell type annotation:
# 1. Finding clusters at different resolutions
# 2. Identifying marker genes for each cluster
# 3. Annotating cell types based on known markers
# 4. Visualizing cell type distribution

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(42)

# Load integrated data
cat("Loading integrated data...\n")
pd <- readRDS("data/integrated_seurat_object.rds")

# Find clusters at different resolutions
cat("Finding clusters at different resolutions...\n")
pd <- FindClusters(object = pd, graph.name = "harmony_snn", resolution = c(0.1, 0.2, 0.3, 0.5))

# Visualize clusters at resolution 0.1
cat("Visualizing clusters...\n")
p1 <- DimPlot(pd, reduction = "umap", group.by = "harmony_snn_res.0.1", label = TRUE)
ggsave("figures/umap_clusters_res0.1.pdf", plot = p1, width = 10, height = 8)

# Find markers for each cluster
cat("Finding marker genes for each cluster...\n")
# Switch to RNA assay for DE testing
DefaultAssay(pd) <- "RNA"
pd <- NormalizeData(pd)
pd.markers <- FindAllMarkers(pd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers to file
write.csv(pd.markers, file = "results/cluster_markers.csv", row.names = FALSE)

# Get top markers for each cluster
top10 <- pd.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create heatmap of top markers
cat("Creating marker gene heatmap...\n")
pd <- ScaleData(pd, features = unique(top10$gene))
p2 <- DoHeatmap(pd, features = top10$gene, group.by = "harmony_snn_res.0.1")
ggsave("figures/marker_heatmap.pdf", plot = p2, width = 12, height = 15)

# Known cell type markers
cat("Visualizing known cell type markers...\n")
cell_markers <- c(
  "GFAP",    # Astrocytes
  "CD74",    # Microglia
  "MOG",     # Oligodendrocytes
  "PDGFRA",  # OPC
  "SLC17A7", # Excitatory neurons
  "GAD1",    # Inhibitory neurons
  "VWF",     # Vascular
  "TRBC2"    # T-cells
)

# Feature plots for known markers
feature_plots <- FeaturePlot(
  pd, 
  features = cell_markers, 
  ncol = 2, 
  reduction = "umap"
)
ggsave("figures/cell_type_markers.pdf", plot = feature_plots, width = 12, height = 16)

# Dot plot for known markers
dot_plot <- DotPlot(pd, features = cell_markers, group.by = "harmony_snn_res.0.1") + 
  RotatedAxis() + 
  ggtitle("Expression of canonical markers by cluster")
ggsave("figures/cell_type_dotplot.pdf", plot = dot_plot, width = 10, height = 8)

# Cell type annotation
cat("Annotating cell types...\n")
# Set cluster identities
Idents(pd) <- "harmony_snn_res.0.1"

# Broad cell type annotation
broad_cell_types <- c(
  "ODC",   # Oligodendrocytes
  "ExN",   # Excitatory Neurons
  "Astro", # Astrocytes
  "ExN",   # Excitatory Neurons
  "MG",    # Microglia
  "InN",   # Inhibitory Neurons
  "OPC",   # Oligodendrocyte Precursor Cells
  "InN",   # Inhibitory Neurons
  "ExN",   # Excitatory Neurons
  "ExN",   # Excitatory Neurons
  "ExN",   # Excitatory Neurons
  "ExN"    # Excitatory Neurons
)
names(broad_cell_types) <- levels(pd)
pd <- RenameIdents(pd, broad_cell_types)
pd@meta.data$BroadCellType <- as.character(Idents(pd))

# Detailed cell type annotation
Idents(pd) <- "harmony_snn_res.0.1"
detailed_cell_types <- c(
  "ODC",   # Oligodendrocytes
  "ExN1",  # Excitatory Neurons 1
  "Astro", # Astrocytes
  "ExN2",  # Excitatory Neurons 2
  "MG",    # Microglia
  "InN1",  # Inhibitory Neurons 1
  "OPC",   # Oligodendrocyte Precursor Cells
  "InN2",  # Inhibitory Neurons 2
  "ExN3",  # Excitatory Neurons 3
  "ExN4",  # Excitatory Neurons 4
  "ExN5",  # Excitatory Neurons 5
  "ExN6"   # Excitatory Neurons 6
)
names(detailed_cell_types) <- levels(pd)
pd <- RenameIdents(pd, detailed_cell_types)
pd@meta.data$CellType <- as.character(Idents(pd))

# Add vascular cells from higher resolution clustering
cat("Adding vascular cell annotation from higher resolution clustering...\n")
pd@meta.data[(pd@meta.data$harmony_snn_res.0.2 == 15), "CellType"] <- "Vas"
pd@meta.data[(pd@meta.data$harmony_snn_res.0.2 == 15), "BroadCellType"] <- "Vas"

# Generate cell type plots
cat("Generating cell type plots...\n")
# Set color scheme
cell_type_colors <- c(
  'Astro' = 'mediumpurple',
  'ExN1' = 'tomato2',
  'ExN2' = 'lightsalmon',
  'ExN3' = 'darkorange',
  'ExN4' = 'gold',
  'ExN5' = 'red',
  'ExN6' = 'brown',
  'InN1' = 'lightblue4',
  'InN2' = 'lightblue3',
  'MG' = 'maroon',
  'ODC' = 'olivedrab',
  'OPC' = 'green',
  'Vas' = 'mistyrose3'
)

# UMAP plots with cell type annotations
Idents(pd) <- "CellType"
p3 <- DimPlot(pd, reduction = "umap", label = TRUE, pt.size = 0.5, cols = cell_type_colors) + 
  ggtitle("Cell Types") + 
  theme(legend.position = "right")
ggsave("figures/umap_cell_types.pdf", plot = p3, width = 10, height = 8)

Idents(pd) <- "BroadCellType"
p4 <- DimPlot(pd, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Broad Cell Types") + 
  theme(legend.position = "right")
ggsave("figures/umap_broad_cell_types.pdf", plot = p4, width = 10, height = 8)

# Cell type distribution by condition
p5 <- ggplot(pd@meta.data, aes(x = SamplID, fill = CellType)) + 
  geom_bar(position = "fill") + 
  RotatedAxis() + 
  theme_classic() + 
  scale_fill_manual(values = cell_type_colors) +
  ggtitle("Cell Type Distribution by Sample") +
  ylab("Proportion") +
  xlab("Sample")
ggsave("figures/cell_type_distribution_by_sample.pdf", plot = p5, width = 12, height = 8)

p6 <- ggplot(pd@meta.data, aes(x = Mutation, fill = CellType)) + 
  geom_bar(position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = cell_type_colors) +
  ggtitle("Cell Type Distribution by Mutation") +
  ylab("Proportion") +
  xlab("Mutation")
ggsave("figures/cell_type_distribution_by_mutation.pdf", plot = p6, width = 10, height = 8)

# Save annotated object
cat("Saving annotated Seurat object...\n")
saveRDS(pd, file = "data/annotated_seurat_object.rds")

cat("Clustering and cell type annotation completed!\n")