# spvAPA v0.1.0 (released on 2024/11/28)
spvAPA: multimodal imputation and supervised analysis of alternative aolyadenylation on single-cell and spatial transcriptome data  

## About  
Alternative polyadenylation (APA) is an essential post-transcriptional modification during messenger RNA (mRNA) maturation in eukaryotic cells. By calculating the relative usage rate of polyA sites, an APA matrix Φ, different from gene expression, can be obtained to interpret gene specificity from a novel perspective. With the development of single-cell and spatial transcriptome sequencing technologies, we have discovered the potential to analyze gene changes at a higher resolution. However, due to the sparsity of APA matrices and the unique nature of the data distribution, there is an urgent need to develop analytical framework suitable for APA matrices. We propose spvAPA for imputation and supervised analysis of APA matrices. spvAPA integrates information from both the gene expression matrix and the Φ matrix to impute the Φ matrix and restore APA signals. Additionally, spvAPA utilizes known prior labels to identify APA features and key APA genes related to the labels through supervised analysis, thereby improving the effectiveness of dimensionality reduction and visualization.  

* The spvAPA package mainly consists of three modules.

<img src="img/Overview.png" width="50%" />  

A. APA usage calculation: Use Poly(A) sites identification tool `scAPAtrap` to extract Poly(A) sites and calculate the usage rate matrix Φ.  

B. APA signal recovery: Integrate multi-modal data from the gene expression and APA usage based on WNN (Weighted Nearest Neighbors) to accurately identify the nearest neighbors for imputing missing entries in the ∅ matrix.  

C. Supervised analysis of APA: Based on sPLS-DA (sparse Partial Least Squares Discriminant Analysis), use prior labels such as histological information to perform supervised analysis of APA information in single-cell and spatial transcriptomics, identifying APA features related to the prior categories and key feature genes. Additionally, spvAPA integrates a flexible visualization component that considers both the selected features and the dual modalities of gene expression and APA, thereby enhancing the visualization of scRNA-seq or spatial transcriptomics data.  

## Getting started  
### Mandatory  
* R (>=4.2.5) (https://www.r-project.org/) is recommended.

### Required R Packages  
* aricode, cluster, clusterProfiler, ClusterR, cowplot, dplyr, ggplot2, magrittr, mclust, mixOmics, RColorBrewer, Seurat

### Installation  
* Install the R package using the following commands on the R console:

```
install.packages("devtools")
require(devtools)
install_github("BMILAB/spvAPA")
library(spvAPA)
browseVignettes('spvAPA')

##or you can download ZIP, and then unzip
devtools::install_local("your_path_of_spvAPA-master.zip", build_vignettes = TRUE)
```

## Application examples  
Vignettes can be found [here](...). Or you can also browse the vignette using the following command on the R console:
```
browseVignettes('spvAPA')
```

### 1. Overview  

An analysis workflow using gene expression matrix and RUD matrix as input, using sPLS-DA and UMAP to perform supervised dimensionality reduction and multimodal integration of scRNA-seq data and spatial transcriptome data.

### 2. Demo data  

Demo data for these vignettes can be downloaded [here](https://github.com/BMILAB/spvAPA/data).  
In the spvAPA package, we provide spatial transcriptomics data from mouse olfactory bulb sample 11 as an example dataset. It includes 260 spatial spots, which are assigned to 5 morphological layers.

```
data("ST11GEM",package = "spvAPA")
data("ST11RUD",package = "spvAPA")
data("ST11label",package = "spvAPA")
dim(ST11GEM)
#> [1] 16218   260
dim(ST11RUD)
#> [1] 4845  260
```

### 3. The imputation of the RUD matrix

WNNimpute() identifies the nearest neighbor cells based on the RUD and gene expression matrices, and then performs iterative imputation using mean filling.

```
RUD_imp <- WNNImpute(gene = ST11GEM,APA = ST11RUD,k = 15,init = F,is.weight = F)
```
### 4. sPLS-DA Supervised Analysis

#### 4.1 sPLS-DA on RUD matrix
The sPLS-DA calculation is implemented through the mixOmics package, which can directly obtain the optimal number of principal components and the number of genes constituting the principal components; the tutorial of model construction can be found in Section 5 of the mixOmics website.

Input: (1) A matrix with rows representing cell/spot barcodes and columns representing genes or other factors; (2) Labels corresponding to cell types or morphological layer for each barcode.

```
RUD_imp_t <- t(RUD_imp)
layer <- ST11label$label
list.keepX <- c(1:10,  seq(20, 150, 10))
RUD_spls <- sPLSDA(X = RUD_imp_t,Y = layer,ncomp = 4,KeepX = list.keepX,folds = 10,invisible = F,nCore = 12)
```
#### 4.2 sPLS-DA on gene expression matrix  

```
GEM_t <- t(ST11GEM)
layer <- ST11label$label
list.keepX <- c(1:10,  seq(20, 150, 10))
GEM_spls <- sPLSDA(X = GEM_t,Y = layer,ncomp = 4,KeepX = list.keepX,folds = 10,invisible = F,nCore = 12)
```

#### 4.3 Select stable genes (variable selection)

During the repeated cross-validation process in perf() we can record how often the same variables are selected across the folds. We select genes with a selection frequency greater than 0.9 as stable feature genes.

```
library(org.Mm.eg.db)
RUD_perf <- RUD_spls[[3]]
featureSelect(perf = RUD_perf,stable = 0.9)
```

### 5. visualisation  

#### 5.1 Spatial distribution of component scores  

Each component could be considered as a meta-gene, and the overall score of each meta-gene is a linear combination of the APA features of the corresponding component. We demonstrated the representation of the principal components in the original spatial coordinates.

```
# Obtain the score matrix of components.
RUD_mod <- RUD_spls[[2]]
comp_score <- as.data.frame(RUD_mod$variates$X)

# Plot

library(cowplot)

comp_score_plot <- PlotComppos(sPLS = RUD_mod,label = ST11label,comps = "ALL",size = 1.5)
plot_grid(plotlist = comp_score_plot, align = "hv")
```

#### 5.2 2D Visualization

##### 5.2.1 Supervised Dimensionality Reduction Visualization for Single Modality

sPLS-DA generates multiple principal components, each containing additional information to distinguish the categories of interest. To capture data patterns arising from multiple known sources of variance, we assume that multiple principal components can be used as inputs for unsupervised dimensionality reduction methods.We performed supervised dimensionality reduction on the gene expression matrix and the RUD matrix using sPLS-DA, and used the principal components from sPLS-DA as feature inputs for UMAP to generate sPLS-DA-UMAP embeddings that effectively separate the categories.

###### RUD  

```
PlotCompumap(spls = RUD_mod,size = 1.7)
```
###### Gene expression matrix

```
GEM_mod <- GEM_spls[[2]]
PlotCompumap(spls = RUD_mod,size = 1.7)
```

##### 5.2.2 Supervised Dimensionality Reduction Visualization Integrating Multimodal Information  

Given that the sPLS-DA-UMAP embedding can separate categories within each label, we tested whether a combined approach for dimensionality reduction could be utilized to visualize gene expression and RUD in a single biplot. We combined the sPLS-DA-UMAP results from two separate sPLS-DA analyses of gene expression and RUD into a single table and used this as the feature set input for UMAP, generating a combined sPLS-DA-UMAP embedding. The resulting biplot preserved biologically meaningful patterns and achieved better intra-class cohesion and inter-class separation.

Compared to the single-modality visualization of the gene expression matrix and the RUD matrix, the two-dimensional visualization combining principal components from both modalities achieved better separation of class labels, successfully distinguishing the GL layer from the OPL layer.

```
PlotCompumap(spls = list(RUD_mod,GEM_mod),size = 1.7)
```





