# Parkinson's Disease Single-Nuclei RNA-Seq Analysis

A comprehensive pipeline for analyzing single-nuclei RNA-seq data from brain samples of Parkinson's disease patients with LRRK2 and GBA mutations compared to healthy controls.

## Overview

This repository contains a collection of R scripts for analyzing single-nuclei RNA-seq data from brain tissue samples of Parkinson's disease patients with different genetic backgrounds (LRRK2, GBA) and healthy controls (HC). The primary analysis includes quality control, normalization, integration and clustering

## Data

The analysis is performed on single-nuclei RNA-seq data from:
- Brain samples (PFC: Prefrontal Cortex and ACC: Anterior Cingulate Cortex)
- Genetic backgrounds: LRRK2 mutation, GBA mutation, and healthy controls
- 12 samples in total

## Pipeline

1. **Data Processing**: 
   - Reading 10X Genomics data
   - Creating Seurat objects
   - Doublet removal using DoubletFinder
   - Quality control filtering

2. **Normalization and Integration**:
   - SCTransform normalization
   - Harmony integration to correct for batch effects
   - PCA and UMAP dimensionality reduction

3. **Clustering and Cell Type Annotation**:
   - Unsupervised clustering
   - Cell type annotation using marker genes
   - Visualization of cell types across different conditions


## Requirements

- R ≥ 4.0
- Seurat ≥ 4.0
- Harmony
- DoubletFinder
- data.table
- ggplot2
- dplyr


## Usage

The repository is organized into modular scripts that can be run sequentially:

1. `01_preprocessing.R`: Data loading and quality control
2. `02_integration.R`: Normalization and integration
3. `03_clustering.R`: Clustering and cell type annotation


