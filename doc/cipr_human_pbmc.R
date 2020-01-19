## ---- include=FALSE-----------------------------------------------------------

knitr::opts_chunk$set(eval=T)

knitr::opts_chunk$set(collapse = TRUE, comment = "#>")



## ---- eval=F------------------------------------------------------------------
#  
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  
#  
#  # Use this option if you want to build vignettes during installation
#  # This can take a long time due to the installation of suggested packages.
#  devtools::install_github(..., build = TRUE, build_opts = c("--no-resave-data", "--no-manual")
#  
#  # Use this if you would like to install the package without vignettes
#  # devtools::install_github("atakanekiz/CIPR-Package")
#  

## ---- eval=T------------------------------------------------------------------

library(dplyr)
library(Seurat)
library(CIPR)


## -----------------------------------------------------------------------------

# Download data

temp <- tempfile()
tempd <- tempdir()

download.file("https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", destfile = temp)

untar(temp, exdir = tempd)

unlink(temp)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = paste0(tempd, "\\filtered_gene_bc_matrices\\hg19"))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc



## -----------------------------------------------------------------------------

# Calculate mitochondrial gene representation (indicative of low quality cells)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter out genes with feature counts outside of 200-2500 range, and >5% mt genes 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## ---- results="hide", message=F-----------------------------------------------

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)


## ---- results="hide", message=F-----------------------------------------------

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


## ---- results="hide", message=F-----------------------------------------------

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


## ---- eval=T------------------------------------------------------------------

ElbowPlot(pbmc)


## ---- results="hide", message=F-----------------------------------------------

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


## ---- results="hide", message=F-----------------------------------------------

pbmc <- RunUMAP(pbmc, dims = 1:10)

pbmc$unnamed_clusters <- Idents(pbmc)


## -----------------------------------------------------------------------------

# saveRDS(pbmc, "pbmc.rds")


## ---- echo=F, results="hide"--------------------------------------------------

allmarkers <- FindAllMarkers(pbmc)


## ---- include=F---------------------------------------------------------------

# saveRDS(allmarkers, "allmarkers.rds")


## ---- results="hide"----------------------------------------------------------

avgexp <- AverageExpression(pbmc)

avgexp <- avgexp$RNA

avgexp$gene <- rownames(avgexp)


## ---- include=F---------------------------------------------------------------

# saveRDS(avgexp, "avgexp.rds")


## ---- include=F---------------------------------------------------------------

# pbmc <- readRDS("pbmc.rds")


## -----------------------------------------------------------------------------

DimPlot(pbmc)


## ---- eval=T, include=F-------------------------------------------------------

# allmarkers <- readRDS("allmarkers.rds")


## ---- eval=T, fig.width=16, fig.height=32, message=F--------------------------

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     # axis.text.x=element_text(color="red") # arguments to pass to ggplot2::theme() to change plotting parameters
     )




## ---- eval=T, fig.width=16, fig.height=5, message=F---------------------------

library(ggplot2)

ind_clu_plots$cluster6 +
    theme(axis.text.y = element_text(color="red"),
          axis.text.x = element_text(color="blue")) +
    labs(fill="Reference")+
    ggtitle("Figure S4a. Automated cluster annotation results are shown for cluster 6") +
    annotate("text", label="2 sd range", x=10, y= 500, size=8, color = "steelblue")+
    annotate("text", label= "1 sd range", x=10, y=175, size=8, color ="orange2")+
  geom_rect(aes(xmin=94, xmax=99, ymin=550, ymax=900), fill=NA, size=3, color="red")
    
    


## ---- eval=T, fig.width=8, fig.height=4.5, message=F--------------------------

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = F,
     plot_top = T, 
     global_results_obj = T, 
     global_plot_obj = T)



## ---- eval=T------------------------------------------------------------------

DT::datatable(CIPR_top_results)


DT::datatable(head(CIPR_all_results))


## ---- eval=T, fig.width=16, fig.height=32, message=F--------------------------

CIPR(input_dat = avgexp,
     comp_method = "all_genes_spearman", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T)




## ---- eval=T, fig.width=8, fig.height=4.5, message=F--------------------------

CIPR(input_dat = avgexp,
     comp_method = "all_genes_spearman", 
     reference = "hsrnaseq", 
     plot_ind = F,
     plot_top = T, 
     global_results_obj = T, 
     global_plot_obj = T)



## ---- eval=T------------------------------------------------------------------

DT::datatable(CIPR_top_results)


DT::datatable(head(CIPR_all_results))


## ---- eval=T, fig.width=16, fig.height=32, message=F--------------------------

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     select_ref_subsets = c("CD4+ T cell", "CD8+ T cell", "Monocyte", "NK cell"))



## ---- eval=T, fig.width=16, fig.height=32, message=F--------------------------

CIPR(input_dat = avgexp,
     comp_method = "all_genes_spearman", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     keep_top_var = 10)




## -----------------------------------------------------------------------------

sessionInfo()



