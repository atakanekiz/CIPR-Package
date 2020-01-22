# CIPR-Package

<br>

<img align="right" src="https://github.com/atakanekiz/CIPR-Package/raw/master/doc/CIPR_hex_mid.png" width=300> 

## Cluster Identity Predictor 

<br>

During the analysis of single cell RNA sequencing (scRNAseq) data, annotating the biological identity of cell clusters is an important step before downstream analyses and it remains technically challenging. The current solutions for annotating single cell clusters generally lack a graphical user interface, can be computationally intensive or have a limited scope. On the other hand, manually annotating single cell clusters by examining the expression of marker genes can be subjective and labor-intensive.

To improve the quality and efficiency of annotating cell clusters in scRNAseq data, we present a web-based R/Shiny app and R package, __Cluster Identity PRedictor (CIPR)__, which provides a graphical user interface to quickly score gene expression profiles of unknown cell clusters against mouse or human references, or a custom dataset provided by the user. CIPR can be easily integrated into the current pipelines to facilitate scRNAseq data analysis.

CIPR performs analyses at individual cluster level and generates informative graphical outputs to help the users assess the quality of algorithmic predictions (see the example outputs below).  


This repository contains the source code for the R package implementation of CIPR pipeline. For CIPR-Shiny, please check out [CIPR-Shiny repository](https://github.com/atakanekiz/CIPR-Shiny).

---

## Installation and Usage

```{r}

install_github("atakanekiz/CIPR-Package", build_vignettes = TRUE)

# # For faster installation without vignette
# install_github("atakanekiz/CIPR-Package", build_vignettes = FALSE)
```

#### Example use case in conjunction with Seurat pipeline

```{r}


library(Seurat)

allmarkers <- FindAllMarkers(seurat_object)
avgexp <- AverageExpression(seurat_object)


# Plot summarizing top scoring references per cluster (logFC comparison)
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T)
     
# Plot summarizing top scoring references per cluster (all-genes correlation)
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T)
     
     
# Plots for individual clusters
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = T,
     plot_top = F)

# Limiting the analysis to certain reference subsets
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T, 
     select_ref_subsets = c("T cell", "B cell", "NK cell"))




```


## Reference datasets available in CIPR

* [Immunological Genome Project (ImmGen)](https://www.immgen.org) microarray data from sorted mouse immune cells. This dataset is prepared by using both V1 and V2 ImmGen releases and it contains 296 samples from 20 different cell types (253 subtypes).

* Mouse RNAseq data from sorted cells reported in [Benayoun et al. (2019)](http://www.genome.org/cgi/doi/10.1101/gr.240093.118). This dataset contains 358 sorted immune and nonimmune samples from 18 different lineages (28 subtypes).

* [Blueprint](https://doi.org/10.3324/haematol.2013.094243)/[Encode](https://doi.org/10.1038/nature11247) RNAseq data that contains 259 sorted human immune and nonimmune samples from 24 different lineages (43 subtypes).

* [Human Primary Cell Atlas](https://doi.org/10.1186/1471-2164-14-632) that contains microarray data from 713 sorted immune and nonimmune cells (37 main cell types and 157 subtypes).

* [DICE (Database for Immune Cell Expression(/eQTLs/Epigenomics)](https://doi.org/10.1016/j.cell.2018.10.022) that contains 1561 human immune samples from 5 main cell types (15 subtypes). To reduce object sizes, mean TPM values per cell type is used.

* Human microarray data from sorted hematopoietic cells reported in [Novershtern et al. (2011)](https://doi.org/10.1016/j.cell.2011.01.004). This dataset contains data from 211 samples and 17 main cell types (38 subtypes)

* Human RNAseq data from sorted cells reported in [Monaco et al. (2019)](https://doi.org/10.1016/j.celrep.2019.01.041). This dataset contains 114 samples originating from 11 main cell types (29 subtypes)

* ___A custom reference dataset provided by the user.___ This dataset can be obtained from a number of high througput methods including microarray and bulk RNAseq. For details about how to prepare custom reference, please see the How-to tab on the [Shiny website](https://aekiz.shinyapps.io/CIPR).

---

## Analytical approach

CIPR calculates pairwise identity scores between individual unknown clusters and the reference samples and generates a vector of identity scores per each cluster in the experiment. While doing this CIPR utilizes two main approaches:

* ___Comparison of differentially expressed genes.___ In this method users provide an input data frame that contains the log fold-change (logFC) values of differentially expressed genes in each cluster. The algorithm first calculates differential expression within the reference data frame for each gene by taking the ratio of the expression value of individual subsets to the average expression in the entire data frame. Then the CIPR pipeline compares these reference logFC values to the logFC from the experimental clusters. The users can select one of three methods for these comparisons:

   * __LogFC dot product:__ LogFC values of the matching genes are mutliplied and added together to yield an aggregate identity score.
   * __LogFC Spearman's correlation:__ Rank correlation is calculated between the logFC values of the experimental and reference data. 
   * __LogFC Pearson's correlation__: Linear correlation is calculated between the logFC values of the expermental and reference data.
   
   <br>
   
* ___Comparison of all genes.___ In this method, users provide an input data frame that contains average gene expression per cluster. The algorithm compares the expression profiles of individual cluster to that from reference dataset. In this method, all the common genes between experimenal and reference data are used in the analysis regardless of their expression values and differential expression status. Users can use one of the two methods in this approach:

   * __Spearman's correlation:__ Rank correlation between the experimental clusters and reference cell subsets
   * __Pearson's correlation:__ Linear correlation (which could be beneficial especially when using custom references where the reference and the experimental data is obtained using similar methodologies.)


---



## Flexible options

To be adaptable to various experimental contexts, CIPR enables users to:

* Select only interesting reference subsets from the provided reference datasets

* Limit the analysis to the genes whose expression variance (in the reference dataset) is above a certain quantile determined by the user.


---

## Sample outputs

### Results per cluster

In the plot below x-axis signifies the individual samples within the reference data frame (ImmGen in this example). Reference cell types are marked by different colors. Each data point indicates the identity score calculated for Cluster 1 in the input data. Shaded regions demarcate 1 and 2 standard deviations around the average identity score across the reference dataset. In this analysis logFC dot product method was used.


<kbd>
<img src=https://github.com/atakanekiz/CIPR-Package/raw/master/doc/sample_ind_output.png>
</kbd>

### Summary of top hits per cluster

It is often easier to examine the top predictions in one graph. This plot shows the top 5 scoring reference samples for each cluster (shown in different colors).

<kbd>
<img src=https://github.com/atakanekiz/CIPR-Package/raw/master/doc/sample_top_output.png>
</kbd>



