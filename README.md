# Parkinson's Disease Single-Nuclei RNA-Seq Analysis

A comprehensive pipeline for analyzing single-nuclei RNA-seq data from brain samples of Parkinson's disease patients with LRRK2 and GBA mutations compared to healthy controls.

## Overview

This repository contains a collection of R scripts for analyzing single-nuclei RNA-seq data from brain tissue samples of Parkinson's disease patients with different genetic backgrounds (LRRK2, GBA) and healthy controls (HC). The primary analysis includes quality control, normalization, integration, and clustering.

## Data

The analysis is performed on single-nuclei RNA-seq data from:
- Brain samples: **PFC (Prefrontal Cortex)** and **ACC (Anterior Cingulate Cortex)**
- Genetic backgrounds: **LRRK2 mutation, GBA mutation, and healthy controls**
- **12 samples** in total

## Pipeline

### 1. Data Processing  
   - Reading **10X Genomics** data  
   - Creating **Seurat** objects  
   - **Doublet removal** using *DoubletFinder*  
   - **Quality control filtering**  

### 2. Normalization and Integration  
   - **SCTransform** normalization  
   - **Harmony** integration to correct for batch effects  
   - **PCA** and **UMAP** dimensionality reduction  

### 3. Clustering and Cell Type Annotation  
   - **Unsupervised clustering**  
   - **Cell type annotation** using **marker genes**  
   - **Visualization** of cell types across different conditions  

## Requirements

Ensure you have the following dependencies installed:  

- **R ≥ 4.0**  
- **Seurat ≥ 4.0**  
- **Harmony**  
- **DoubletFinder**  
- **data.table**  
- **ggplot2**  
- **dplyr**  

To install the required R packages, run:  
```r
install.packages(c("Seurat", "data.table", "ggplot2", "dplyr"))
remotes::install_github("immunogenomics/harmony")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

```

## Data Availability

The single-nuclei RNA-seq data used in this analysis is publicly available in the **Gene Expression Omnibus (GEO)** under the accession number:  

**GSE272760**  

### Downloading the Data  

You can download the data using the `GEOquery` R package:  

```r
# Install GEOquery if not installed
if (!requireNamespace("GEOquery", quietly = TRUE)) install.packages("GEOquery")

# Load the package
library(GEOquery)

# Download the dataset
gse <- getGEO("GSE272760", GSEMatrix = TRUE, getGPL = FALSE)

# View sample metadata
metadata <- pData(gse[[1]])
head(metadata)




