# spvAPA v0.1.0 (released on 2024/11/28)
Supervised analysis of alternative polyadenylation from single-cell and spatial transcriptomics data with spvAPA

## About  
spvAPA is a supervised analytical framework specifically designed for alternative polyadenylation (APA) analysis from both single-cell and spatial transcriptomics data. Firstly, an iterative imputation method based on weighted nearest neighbor (WNN) was designed for imputing missing entries in the APA data (∅ matrix). Secondly, a supervised feature selection method based on sparse partial least squares discriminant analysis (sPLS-DA) is devised to identify APA features from scRNA-seq and spatial transcriptomics data. Additionally, spvAPA integrates a flexible visualization module that considers both the selected features and the dual modalities of gene expression and APA, thereby enhancing the visualization of high dimensional scRNA-seq or spatial transcriptomics data.

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
Vignettes can be found [here](https://github.com/BMILAB/spvAPA/blob/master/doc/spvAPA_vignette.html). Or you can also browse the vignette using the following command on the R console:
```
browseVignettes('spvAPA')
```
### Identification of APA sites and calculation of the APA ∅ matrix  

scAPAtrap was used to identify APA sites from single-cell RNA-seq and spatial transcriptomics data. Then movAPA was used to calculate the APA ∅ matrix based on the RUD (relative usage of distal poly(A) site) metric.
The following provides example code for users to refer to on how to obtain the APA matrix from their own data.

[scAPAtrap](https://github.com/BMILAB/scAPAtrap) was used to identify APA sites from single-cell RNA-seq and spatial transcriptomics data. Then movAPA was used to calculate the APA ∅ matrix based on the RUD (relative usage of distal poly(A) site) metric.
The following provides example code for users to refer to on how to obtain the APA matrix from their own data.

#### Get PA count matrix  from BAM file with the scAPAtrap  

```R
library(scAPAtrap)

## input BAM
dir0='/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/xwu/demoFly/'
inputBam=paste0(dir0, "fly_demo.bam")
logf=gsub('.bam', '.APA.notails.onestep.log', inputBam, fixed=TRUE)
outputDir="APA.notails.onestep"
## full path of tools
tools=list(samtools='/home/dell/miniconda3/envs/scAPA/bin/samtools',
umitools='/home/dell/miniconda3/envs/scAPA/bin/umi_tools',
featureCounts="/home/dell/miniconda3/envs/scAPA/bin/featureCounts",
star='/home/dell/miniconda3/envs/scAPA/bin/STAR')
## set parameters, first load default parameters and then modify
trap.params=setTrapParams()
trap.params$tails.search='no'
trap.params$chrs=c('2L','2R','3L','3R','4','X','Y')
## barcode (if have)
trap.params$barcode <- read.delim2(paste0(dir0, 'barcode.txt'), header = F)$V1
## Run scAPAtrap
scAPAtrap(tools=tools,
trap.params=trap.params,
inputBam=inputBam,
outputDir=outputDir,
logf=logf)
```



#### Get APA ∅ matrix from scAPAtrap results with the movAPA package



This dataset contains many PAs (also called peaks in scAPAtrap), which may not be suitable for downstream analysis. First, we can create a PACdataset object and remove extremely lowly expressed peaks. Here we suggest removing PAs with <10 tags supported by all cells and PAs that are expressed in <10 cells.

```R
library(movAPA)
PACds=createPACdataset(counts=scAPAtrapData$peaks.count, anno=scAPAtrapData$peaks.meta) # create a PACdataset object.
PACds=subsetPACds(PACds, totPACtag = 10, minExprConds = 10, verbose=TRUE)
```

After read the data into a PACdataset, users can use many functions in movAPA for removing internal priming artifacts, annotating PACs, polyA signal analysis, etc. The process can be referred to in the movAPA package.

##### Remove internal priming artifacts  

```R
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE, verbose = FALSE)
bsgenome<-BSgenome.Hsapiens.UCSC.hg38
PACdsIP=removePACdsIP(PACds, bsgenome, returnBoth=TRUE,
up=-10, dn=10, conA=6, sepA=7)
PACds=PACdsIP$real
```

##### Annotate genomic regions for PACs  

Genes with or without annotated 3’UTR could be assigned an extended 3’UTR of a given length using the function ext3UTRPACds, which can improve the “recovery” of poly(A) sites falling within authentic 3’UTRs. For single cell data, we suggest analyzing only 3’UTR PAs and discarding PAs from other regions.

```R
library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = TRUE, verbose = FALSE)
txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
# annotate PAs
PACds=annotatePAC(PACds, txdb)
# extend 3UTR by 1000bp
PACds=ext3UTRPACds(PACds, 1000)
```

##### Get APA ∅ matrix using the smart RUD method  

The smartRUD indicator is provided in movAPA v0.2.0, and it is recommended to use it. Pay attention to setting clearPAT=1 to remove cases of PAs with only 1 read count. At the same time, check the distance between the two PAs first to select a suitable dist for filtering the proximal and distal PAs of 3’UTR.

```R
pd=get3UTRAPApd(pacds=PACds, minDist=50, maxDist=5000, minRatio=0.05,
fixDistal=FALSE, addCols='pd')
rud=movAPAindex(pd, method="smartRUD", sRUD.oweight=TRUE, clearPAT=1)
```

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

<figure>
  <img src="img/The spatial representation of principal components.png" width="50%"> 
  <figcaption>The spatial representation of principal components.</figcaption>
</figure>

#### 5.2 2D Visualization

##### 5.2.1 Supervised Dimensionality Reduction Visualization for Single Modality

sPLS-DA generates multiple principal components, each containing additional information to distinguish the categories of interest. To capture data patterns arising from multiple known sources of variance, we assume that multiple principal components can be used as inputs for unsupervised dimensionality reduction methods.We performed supervised dimensionality reduction on the gene expression matrix and the RUD matrix using sPLS-DA, and used the principal components from sPLS-DA as feature inputs for UMAP to generate sPLS-DA-UMAP embeddings that effectively separate the categories.

###### RUD  

```
PlotCompumap(spls = RUD_mod,size = 1.7)
```

<figure>
  <img src="img/2D visualization of RUD single-modality.png" width="50%"> 
  <figcaption>2D visualization of RUD single-modality.</figcaption>
</figure>

###### Gene expression matrix

```
GEM_mod <- GEM_spls[[2]]
PlotCompumap(spls = RUD_mod,size = 1.7)
```
<figure>
  <img src="img/2D visualization of GEM single-modality.png" width="50%"> 
  <figcaption>2D visualization of GEM single-modality.</figcaption>
</figure>

##### 5.2.2 Supervised Dimensionality Reduction Visualization Integrating Multimodal Information  

Given that the sPLS-DA-UMAP embedding can separate categories within each label, we tested whether a combined approach for dimensionality reduction could be utilized to visualize gene expression and RUD in a single biplot. We combined the sPLS-DA-UMAP results from two separate sPLS-DA analyses of gene expression and RUD into a single table and used this as the feature set input for UMAP, generating a combined sPLS-DA-UMAP embedding. The resulting biplot preserved biologically meaningful patterns and achieved better intra-class cohesion and inter-class separation.

Compared to the single-modality visualization of the gene expression matrix and the RUD matrix, the two-dimensional visualization combining principal components from both modalities achieved better separation of class labels, successfully distinguishing the GL layer from the OPL layer.

```
PlotCompumap(spls = list(RUD_mod,GEM_mod),size = 1.7)
```
<figure>
  <img src="img/2D visualization of multi-modality.png" width="50%"> 
  <figcaption>2D visualization of multi-modality.</figcaption>
</figure>




