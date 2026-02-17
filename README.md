# Microbiome Analysis of Colorectal Cancer (YACHIDA_CRC_2019)

This repository contains R scripts for performing a comprehensive microbiome data analysis based on the **YACHIDA_CRC_2019** dataset. The analysis focuses on the taxonomic profile of patients with Colorectal Cancer (CRC) compared to healthy controls, utilizing genus-level count data.

---

## Repository Structure

### 1. `ex1.R`: Data Acquisition & Study Design
This script handles the initial setup and exploratory data analysis (EDA).
* **Data Loading**: Automatically fetches metadata and genus-level count matrices from remote repositories.
* **Demographics**: Calculates the total number of samples and the distribution across clinical study groups, such as "Healthy" and various CRC stages.
* **Visualization**: 
    * Creates a **Barplot** of sample counts per group to visualize the study design.
    * Generates an **Age Density Plot** to compare the age distribution across different health conditions.

### 2. `ex2_and_3.R`: Ecological Diversity & Statistics
This script performs advanced ecological modeling, specifically comparing the **Healthy** group against **Stage_III_IV** patients.

#### **Alpha Diversity (Within-Sample)**
* **Metrics**: Computes **Observed Richness** (total genera count) and the **Shannon Diversity Index** using the `vegan` package.
* **Statistical Testing**: Employs the **Wilcoxon Rank-Sum Test** to identify significant differences in microbial richness and diversity between the two groups.
* **Visualization**: Produces boxplots with integrated p-values for rapid statistical interpretation.

#### **Beta Diversity (Between-Sample)**
* **Dissimilarity**: Calculates **Jaccard Dissimilarity** based on a presence/absence (binary) matrix of the microbial taxa.
* **Ordination**: Performs **Principal Coordinates Analysis (PCoA)** to visualize community clustering.
* **PERMANOVA**: Conducts a Permutational Multivariate Analysis of Variance to determine the effect size ($R^2$) and statistical significance ($p$-value) of the clinical group on the microbiome composition.

---

## Requirements
The analysis relies on the following R packages:
* `tidyverse`: For data manipulation and visualization.
* `vegan`: For ecological diversity indices and PERMANOVA.
* `ggpubr`: For adding statistical layers to plots.
* `patchwork`: For combining multiple plots into single figures.

## How to Run
1. Ensure you have an active internet connection to download the datasets directly from the provided URLs.
2. Run `ex1.R` first to understand the data structure and population demographics.
3. Run `ex2_and_3.R` to execute the statistical comparisons and diversity analyses.
