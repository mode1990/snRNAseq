# 01_preprocessing.R
# 
# This script handles the initial preprocessing steps for single-cell RNA-seq data:
# 1. Reading the 10X Genomics data
# 2. Creating Seurat objects
# 3. Performing initial QC
# 4. Identifying and removing doublets using DoubletFinder
# 5. Merging the datasets

# Load required libraries
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Function to create Seurat object with proper metadata
create_seurat_with_metadata <- function(data_dir, sample_id, mutation, region) {
  # Read 10X data
  data <- Read10X(data.dir = data_dir)
  
  # Create Seurat object with basic filtering
  seurat_obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
  
  # Add metadata
  seurat_obj@meta.data$Mutation <- mutation
  seurat_obj@meta.data$Region <- region
  seurat_obj@meta.data$SamplID <- sample_id
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Basic QC filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1000 & percent.mt < 5)
  
  return(seurat_obj)
}

# Function to remove doublets using DoubletFinder
remove_doublets <- function(seurat_obj) {
  # Pre-process for DoubletFinder
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Get annotations for homotypic doublet estimation
  annotations <- seurat_obj@meta.data$orig.ident
  homotypic.prop <- modelHomotypic(annotations)
  
  # Calculate expected doublet rate (adjust 0.012 based on your 10X chemistry)
  nExp_poi <- round(0.012 * nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = 0.09, 
                                 nExp = nExp_poi, reuse.pANN = FALSE)
  
  # Get the correct column name (it varies based on parameters)
  df_col <- colnames(seurat_obj@meta.data)[grep("DF.classifications", colnames(seurat_obj@meta.data))]
  
  # Subset to keep only singlets
  seurat_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[which(seurat_obj@meta.data[[df_col]] == "Singlet")])
  
  return(seurat_obj)
}

# Read and preprocess all samples
# Sample 1 - GBA mutation in PFC
s1.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample1")
s1 <- CreateSeuratObject(counts = s1.data, min.cells = 3, min.features = 200)
s1@meta.data$Mutation <- "GBA"
s1@meta.data$Region <- "PFC"
s1@meta.data$SamplID <- "Sample1"

# Sample 2 - GBA mutation in PFC
s2.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample2")
s2 <- CreateSeuratObject(counts = s2.data, min.cells = 3, min.features = 200)
s2@meta.data$Mutation <- "GBA"
s2@meta.data$Region <- "PFC"
s2@meta.data$SamplID <- "Sample2"

# Sample 3 - LRRK2 mutation in PFC
s3.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample3")
s3 <- CreateSeuratObject(counts = s3.data, min.cells = 3, min.features = 200)
s3@meta.data$Mutation <- "LRRK2"
s3@meta.data$Region <- "PFC"
s3@meta.data$SamplID <- "Sample3"

# Sample 4 - LRRK2 mutation in PFC
s4.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample4")
s4 <- CreateSeuratObject(counts = s4.data, min.cells = 3, min.features = 200)
s4@meta.data$Mutation <- "LRRK2"
s4@meta.data$Region <- "PFC"
s4@meta.data$SamplID <- "Sample4"

# Sample 5 - GBA mutation in ACC
s5.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample5")
s5 <- CreateSeuratObject(counts = s5.data, min.cells = 3, min.features = 200)
s5@meta.data$Mutation <- "GBA"
s5@meta.data$Region <- "ACC"
s5@meta.data$SamplID <- "Sample5"

# Sample 6 - GBA mutation in ACC
s6.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample6")
s6 <- CreateSeuratObject(counts = s6.data, min.cells = 3, min.features = 200)
s6@meta.data$Mutation <- "GBA"
s6@meta.data$Region <- "ACC"
s6@meta.data$SamplID <- "Sample6"

# Sample 7 - HC (Healthy Control) in PFC
s7.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample7")
s7 <- CreateSeuratObject(counts = s7.data, min.cells = 3, min.features = 200)
s7@meta.data$Mutation <- "HC"
s7@meta.data$Region <- "PFC"
s7@meta.data$SamplID <- "Sample7"

# Sample 8 - HC (Healthy Control) in PFC
s8.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample8")
s8 <- CreateSeuratObject(counts = s8.data, min.cells = 3, min.features = 200)
s8@meta.data$Mutation <- "HC"
s8@meta.data$Region <- "PFC"
s8@meta.data$SamplID <- "Sample8"

# Sample 9 - LRRK2 mutation in ACC
s9.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample9")
s9 <- CreateSeuratObject(counts = s9.data, min.cells = 3, min.features = 200)
s9@meta.data$Mutation <- "LRRK2"
s9@meta.data$Region <- "ACC"
s9@meta.data$SamplID <- "Sample9"

# Sample 10 - LRRK2 mutation in ACC
s10.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample10")
s10 <- CreateSeuratObject(counts = s10.data, min.cells = 3, min.features = 200)
s10@meta.data$Mutation <- "LRRK2"
s10@meta.data$Region <- "ACC"
s10@meta.data$SamplID <- "Sample10"

# Sample 11 - HC (Healthy Control) in ACC
s11.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample11")
s11 <- CreateSeuratObject(counts = s11.data, min.cells = 3, min.features = 200)
s11@meta.data$Mutation <- "HC"
s11@meta.data$Region <- "ACC"
s11@meta.data$SamplID <- "Sample11"

# Sample 12 - HC (Healthy Control) in ACC
s12.data <- Read10X(data.dir = "/path/to/the/filteredmatrix/sample12")
s12 <- CreateSeuratObject(counts = s12.data, min.cells = 3, min.features = 200)
s12@meta.data$Mutation <- "HC"
s12@meta.data$Region <- "ACC"
s12@meta.data$SamplID <- "Sample12"

# Process each sample (example for s1, repeat for all samples)
# Initial QC
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")
s1 <- subset(s1, subset = nFeature_RNA > 1000 & percent.mt < 5)

# Doublet removal (example for s1, repeat for all samples)
# Run DoubletFinder
annotations <- s1@meta.data$orig.ident
homotypic.prop <- modelHomotypic(annotations)        
nExp_poi <- round(0.012*nrow(s1@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
s1 <- doubletFinder_v3(s1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
# Extract the DoubletFinder column name (varies by parameters)
df_cols <- grep("DF.classifications", colnames(s1@meta.data), value = TRUE)
s1 <- subset(s1, cells=rownames(s1@meta.data)[which(s1@meta.data[[df_cols]] == "Singlet")])

# Repeat for all samples s2-s12

# Merge all objects
pd <- merge(s1, 
            y = c(s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12), 
            add.cell.ids = c("GBA1", "GBA2", "LRRK21", "LRRK22", "GBA3", "GBA4", 
                             "HC1", "HC2", "LRRK23", "LRRK24", "HC3", "HC4"), 
            merge.data = TRUE, 
            project = "PD_scRNA")

# Final QC for the merged object
pd[["percent.mt"]] <- PercentageFeatureSet(pd, pattern = "^MT-")
pd[["percent.ribo"]] <- PercentageFeatureSet(pd, pattern = "^RP[SL]")
pd <- subset(pd, subset = nFeature_RNA > 1000 & percent.mt < 5 & percent.ribo < 5)

# QC visualization 
VlnPlot(pd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

plot1 <- FeatureScatter(pd, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_combined <- plot1 + plot2

# Save QC plots
ggsave("figures/qc_violin_plots.pdf", plot = VlnPlot(pd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4), width = 12, height = 6)
ggsave("figures/qc_feature_scatter.pdf", plot = plot_combined, width = 12, height = 6)

# Save the merged Seurat object
saveRDS(pd, file = "/path/to/the/output/merged_seurat_object_preprocessed.rds")

cat("Preprocessing completed and saved to output/merged_seurat_object_preprocessed.rds\n")
