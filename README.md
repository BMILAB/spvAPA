# spvAPA v0.1.0 (released on 2024/05/28)
spvAPA: multimodal imputation and supervised analysis of alternative aolyadenylation on single-cell and spatial transcriptome data  

## About  
Alternative polyadenylation (APA) is an essential post-transcriptional modification during messenger RNA (mRNA) maturation in eukaryotic cells. By calculating the relative usage rate of polyA sites, an APA matrix Φ, different from gene expression, can be obtained to interpret gene specificity from a novel perspective. With the development of single-cell and spatial transcriptome sequencing technologies, we have discovered the potential to analyze gene changes at a higher resolution. However, due to the sparsity of APA matrices and the unique nature of the data distribution, there is an urgent need to develop analytical framework suitable for APA matrices. We propose spvAPA for imputation and supervised analysis of APA matrices. spvAPA integrates information from both the gene expression matrix and the Φ matrix to impute the Φ matrix and restore APA signals. Additionally, spvAPA utilizes known prior labels to identify APA features and key APA genes related to the labels through supervised analysis, thereby improving the effectiveness of dimensionality reduction and visualization.  

* The spvAPA package mainly consists of three modules.

<img src="img/Overview.png" width="60%" />  

