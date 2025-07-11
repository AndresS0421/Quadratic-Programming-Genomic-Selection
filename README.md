# Genomic Selection Optimization Using Quadratic Programming

## Overview

This repository contains the implementation and analysis of a novel genomic selection approach that employs quadratic programming (QP) optimization for optimal breeding line selection. The research compares traditional selection methods with QP-based approaches across multiple crop datasets, focusing on maximizing breeding values while maintaining genetic diversity.

## Research Objectives

- **Primary Goal**: Develop and evaluate a quadratic programming framework for optimal genomic selection
- **Secondary Goals**:
  - Compare QP-based selection with traditional top-performing line selection
  - Evaluate the trade-off between breeding value maximization and genetic diversity maintenance
  - Assess performance across diverse crop species and environments
  - Quantify the impact of different selection intensities and penalty parameters

## Methodology

### Quadratic Programming Formulation

The optimization problem is formulated as:

```
Maximize: t(d_vec) %*% x - k * quad_form(x, Q1)
Subject to:
- x ∈ {0,1}^n (binary selection variables)
- Σx = p*n (select p proportion of total lines)
- x ≥ 0, x ≤ 1
```

Where:
- `d`: vector of BLUEs (Best Linear Unbiased Estimators)
- `D`: genomic relationship matrix (Geno matrix)
- `k`: penalty parameter controlling diversity emphasis
- `p`: proportion of lines to select (0.1 or 0.2)
- `x`: binary selection vector

### Key Features

1. **Multi-objective Optimization**: Balances breeding value maximization with genetic diversity maintenance
2. **Genomic Relationship Integration**: Uses GRM to quantify genetic similarities between lines
3. **Flexible Selection Intensity**: Supports different selection proportions (10% and 20%)
4. **Multiple Similarity Metrics**: Implements cosine similarity, Euclidean distance, Pearson correlation, and Manhattan distance
5. **Comprehensive Evaluation**: Compares selection methods across multiple performance metrics

## Datasets

The research evaluates the following crop datasets:

- **Rice**: Indica and Japonica varieties, IRRI Philippines collection
- **Maize**: Agricultural lines dataset  
- **Eucalyptus**: Australian Calister 2022 collection
- **Groundnut**: Agricultural lines dataset
- **Pinus**: Florida USA Resende 2012 collection

Each dataset contains:
- Phenotypic observations for multiple traits
- Genomic relationship matrices
- Environmental information (when applicable)

## Project Structure

```
QuadraticProgramming/
├── README.md
├── QP_selection.R              # Main optimization implementation
├── utils.R                     # Utility functions
├── Measure_Similarities.R      # Genetic similarity calculations
├── Datasets/                   # Genomic datasets (.RData files)
├── Optimal_Line_Selection/     # Results and output files
├── Results_Graphics/           # Generated plots and visualizations
└── Other_Codes/               # Additional analysis scripts
    ├── Graphicate_Results.R
    ├── Summarise_Selection_Results.R
    ├── MeanRatiosDatasets.R
    └── Filter_K_Values.R
```

## Installation and Dependencies

### Required R Packages

```r
# Genomic analysis
library(BGLR)
library(SKM)

# Data manipulation
library(dplyr)
library(tidyr)
library(reshape2)
library(plyr)

# Optimization
library(quadprog)
library(CVXR)
library(ROI)
library(ROI.plugin.glpk)
library(pracma)

# Visualization
library(ggplot2)
library(caret)
```

### Installation

```r
# Install required packages
install.packages(c("BGLR", "SKM", "dplyr", "tidyr", "reshape2", 
                   "quadprog", "CVXR", "ROI", "ROI.plugin.glpk", 
                   "pracma", "ggplot2", "caret"))
```

## Usage

### Running the Main Analysis

1. **Configure Parameters**: Edit the variables in `QP_selection.R`:
   ```r
   datasets <- c("Indica_AL", "Japonica_AL", "Groundnut_AL", ...)
   dataset_selected_index <- 8  # Select dataset index
   percentages <- c(0.1, 0.2)   # Selection intensities
   k_list <- c(1)               # Penalty parameters
   ```

2. **Execute Main Analysis**:
   ```r
   source("QP_selection.R")
   ```

3. **Generate Visualizations**:
   ```r
   source("Other_Codes/Graphicate_Results.R")
   ```

### Key Functions

- `best_lines_match()`: Calculates percentage matching between selected and top-performing lines
- `cosine_similarity()`: Computes cosine similarity between genomic profiles
- `euclidean_distance()`: Calculates normalized Euclidean distance
- `pearson_correlation()`: Measures linear genetic relationships
- `manhattan_distance()`: Computes Manhattan distance between lines

## Results and Outputs

### Performance Metrics

The analysis generates comprehensive results including:

- **Percentage Matching (PM)**: Overlap between QP-selected and traditionally selected lines
- **Total Breeding Value (TBV)**: Sum of breeding values for selected lines (QP vs Traditional)
- **Mean Breeding Value**: Average breeding value of selected lines (QP vs Traditional)
- **Variance Metrics**: Genetic variance within selected populations (QP vs Traditional)
- **Expected Loss**: Risk-adjusted performance measures (QP vs Traditional)
- **Ratio Metrics**: Gain-to-variance ratios (QP vs Traditional)
- **Average Relatedness**: Genetic similarity within selected groups (QP vs Traditional)

### Output Files

The results are saved in the `Optimal_Line_Selection/[Dataset]/` directory structure:

- `summary_ALL_[dataset].csv`: Comprehensive performance metrics for each dataset
- `selection_ALL_[dataset].csv`: Detailed selection results per line for each dataset
- `summary_GROUP_[dataset].csv`: Grouped summary results across datasets
- Generated plots in `Results_Graphics/[Dataset]/` directories

### Visualization

The project generates various plots comparing:
- QP vs. traditional selection performance
- Ratio comparisons across datasets
- Performance across different selection intensities
- Trade-offs between breeding value and genetic diversity

## Key Findings

1. **Selection Efficiency**: QP optimization can identify lines with comparable or superior breeding values compared to traditional methods
2. **Genetic Diversity**: The penalty parameter `k` effectively controls the diversity-performance trade-off
3. **Dataset Variability**: Performance gains vary across crop species and trait characteristics
4. **Selection Intensity Impact**: Different selection proportions (10% vs 20%) show varying optimization benefits

## Applications

This research framework can be applied to:
- Plant breeding programs seeking optimal parent selection
- Animal breeding for balanced genetic improvement
- Conservation breeding programs requiring diversity maintenance
- Multi-trait selection scenarios
- Cross-population breeding decisions

## Contributing

This research is part of ongoing genomic selection optimization studies. For questions or collaborations, please contact the research team.

## Acknowledgments

We thank the data providers and research institutions that contributed the genomic datasets used in this study. 