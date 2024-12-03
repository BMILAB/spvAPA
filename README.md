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

### Impute missing entries in the APA ∅ matrix  

This [tutorial](https://github.com/BMILAB/spvAPA/blob/master/doc/spvAPA_vignette.html) introduces how to use the gene expression matrix and APA matrix as inputs to obtain the imputed APA ∅ matrix. 
We used single-cell RNA-seq data from mouse olfactory bulb for demonstration. The demo data including the gene expression matrix and APA matrix can be downloaded [here](https://github.com/BMILAB/spvAPA/tree/master/data).

### Application of spvAPA on single-cell RNA-seq data: feature selection, dimensionality reduction, and visualization  

This [tutorial]() introduces how to use sPLS-DA module in spvAPA to identify APA features related to cell types and enhance dimensionality reduction and visualization. 
The input data include the APA matrix after imputation and cell type annotations of the scRNA-seq data from mouse olfactory bulb, which can be downloaded here.

### Application of spvAPA on spatial transcriptomics data: feature selection, dimensionality reduction, and visualization  

This [tutorial](https://github.com/BMILAB/spvAPA/blob/master/doc/spvAPA_vignette.html) introduces how to use sPLS-DA module in spvAPA for feature selection, dimensionality reduction, and visualization on spatial transcriptomics data.
The spatial transcriptomics data of mouse olfactory bulb was used for demonstration. The input data include the APA matrix, gene expression matrix, coordinates of spots, and histological annotations, which can be downloaded [here](https://github.com/BMILAB/spvAPA/tree/master/data).

