#' @title Imputation APA index using gene expression matrix and APA matrix.
#' @description A KNN-based and multimodal imputation model for recovering APA signals in the APA index matrix.
#' @param gene Reference gene expression matrix, with genes as rows and cells or spots as columns.
#' @param APA APA index to be recovered, with genes as rows and cells or spots as columns.
#' @param k Defines k for the k-nearest neighbor algorithm.
#' @param init Whether to initialize the APA index.
#'
#' @return Imputed APA index
#' @export
#' @import magrittr
#' @import Seurat
#'

WNNImpute = function(gene,APA, k = 20, init = TRUE){
  obj.ge.pdui.wnn = CreateSeuratObject(counts = gene)
  APA_noZero <- APA
  APA_noZero[is.na(APA_noZero) ] = 0

  temp<- CreateSeuratObject(counts = APA_noZero)
  obj.ge.pdui.wnn[["PDUI"]]=temp@assays$RNA

  ##Gene Matrix标准化，pca降维，umap
  DefaultAssay(obj.ge.pdui.wnn) <- "RNA"
  obj.ge.pdui.wnn <- SCTransform(obj.ge.pdui.wnn, verbose = T,assay = "RNA",
                                 new.assay.name = "rna_SCT") %>%
    RunPCA(reduction.name = "rna_pca") %>%
    RunUMAP(reduction = "rna_pca",dims = 1:50,
            reduction.name = 'umap.rna',
            reduction.key = 'rnaUMAP_')
  ##RUD标准化
  DefaultAssay(obj.ge.pdui.wnn) <- "PDUI"
  obj.ge.pdui.wnn <-
    FindVariableFeatures(obj.ge.pdui.wnn, selection.method = "vst")%>%
    NormalizeData()%>%
    ScaleData()%>%
    RunPCA(reduction.name = "pa_pca") %>%
    RunUMAP(reduction = "pa_pca",dims = 1:50,
            reduction.name = 'umap.pa',
            reduction.key = 'paUMAP_')
  ## Construct weighted nearest neighbor graph
  obj.ge.pdui.wnn <- FindMultiModalNeighbors(obj.ge.pdui.wnn, reduction.list = list("rna_pca", "pa_pca"), dims.list = list(1:50, 1:50))
  weight=obj.ge.pdui.wnn@neighbors[["weighted.nn"]]@nn.idx[,1:k]
  rownames(weight) <-obj.ge.pdui.wnn@assays[["RNA"]]$counts@Dimnames[[2]]
  colnames(weight) <- paste0("N", 1:k)

  colNames <- colnames(APA)[colnames(APA) %in% colnames(gene)]
  APA <- APA[, colNames]
  print(paste0("NA values before imputed: ", sum(is.na(APA))))
  gene <- gene[, colNames]

  if (init == TRUE) {
    gene_some <- gene[rownames(APA), ]
    sum(gene_some==0)
    is_zero <- gene_some == 0
    APA[is_zero & is.na(APA)] = 0
    print(paste0("NA values after initial imputed: ", sum(is.na(APA))))
  }
  num = 1
  while (sum(is.na(APA)) > 0 & num <= 10) {
    print(paste0("imputing....................", num))
    tmp=APA
    for (i in 1:ncol(tmp)) {
      target <- tmp[, i] # 第i个细胞下每个基因的RUD表达量：
      if(sum(is.na(target)) == 0) next
      tianchong <- tmp[, as.numeric(weight[i, ])] # 第i个细胞最近的k个细胞的RUD表达量:gene x k
      tianchong.mean <- apply(tianchong,1,mean,na.rm = TRUE)
      target[is.na(target)] <- tianchong.mean[is.na(target)]
      APA[, i] <- target
    }
    print(paste0("NA values : ", sum(is.na(APA))))
    num = num + 1
  }
  if (num > 10) {
    print("over the max impute times,convert all Na to Zero")
    APA[is.na(APA)] = 0
  }
  return(APA)
}

#' @title Standard workflow for dimensionality reduction and clustering of the RUD matrix/gene express matrix using the Seurat package.
#' @param seurat A Seurat object.
#' @param norm Whether to normalize the matrix in the assay.
#' @param res Value of the resolution parameter in the \code{FindClusters()} function of the \code{Seurat} package.
#' @param dims Dimensions of reduction to use as input in the \code{FindNeighbors()} function of the \code{Seurat} package.
#' @param nfeatures Number of features to select as top variable features in the \code{FindVariableFeatures()} function of the \code{Seurat} package.
#' @param k.param Defines k for the k-nearest neighbor algorithm in the \code{FindNeighbors()} function of the \code{Seurat} package.
#' @param algorithm Algorithm for modularity optimization in the \code{FindClusters()} function of the \code{Seurat} package.
#' @param normalization.method Method for normalization in the \code{NormalizeData()} function of the \code{Seurat} package.
#'
#' @return A Seurat object that has been reduced dimensionality using PCA and UMAP.
#' @export

makeCluster = function(seurat,norm=T,res=0.1,dims=1:50,nfeatures = 2000,k.param = 30,algorithm = 1,normalization.method='LogNormalize'){
  if(norm){
    seurat <- NormalizeData(seurat, normalization.method = normalization.method,
                            scale.factor = 10000)
  }

  seurat <- FindVariableFeatures(seurat, selection.method = "vst",
                                 nfeatures = nfeatures)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
  seurat <- FindNeighbors(seurat, dims = dims, k.param = k.param)
  seurat <- FindClusters(seurat, resolution = res,algorithm = algorithm)
  seurat <- RunUMAP(seurat,dims = dims)
  return(seurat)
}


#' @title Calculate silhouette coefficient and external metrics (ARI、NMI、Jaccard & Purity) for clustering based on true labels.
#'
#' @param sc A Seurat object.
#' @param true_labels Cell type or morphology labels, etc., must be equal in number to the cells or spots in /{sc}.
#' @import ClusterR
#' @importFrom mclust adjustedRandIndex
#' @import cluster
#' @import aricode
#' @import Seurat
#' @import dplyr
#' @return A list.
#' @export

estimate <- function(sc,true_labels){
  true_labels=as.numeric(as.factor(true_labels))
  reduc_dat = as.data.frame.matrix(Embeddings(sc[['umap']]))
  sil1 = silhouette(true_labels, dist(reduc_dat))
  sil = group_by(as.data.frame.matrix(sil1),cluster)
  sil = as.data.frame(summarise(sil, mean_sil=mean(sil_width)))
  sil = colMeans(sil)[2]
  ################ cluster louvain
  n_class = length(unique(true_labels))
  print(n_class)
  jaccard <- external_validation(true_labels,as.numeric(as.factor(sc$seurat_clusters)), method = "jaccard_index")
  ari_lov = adjustedRandIndex(sc$seurat_clusters,true_labels)
  nmi_lov = NMI(as.vector(sc$seurat_clusters),as.vector(true_labels))
  c2 = ClusteringPurity(sc$seurat_clusters,as.matrix(true_labels))
  return(list(sc=sil,ari=ari_lov,nmi=nmi_lov,purity=c2,jaccard=jaccard))
}


#' @title Run sPLS (sparse Partial Least Squares) supervised analysis.
#'
#' @param X Matrix, rows are cells or spots, columns are genes.
#' @param Y Annotations of cells or spots corresponding to barcodes, a vector.
#' @param ncomp The expected number of components to be found.
#' @param KeepX  The gene composition of each component.Grid of possible keepX values that will be tested for each comp
#' @param folds The number of folds for cross-validation in classification.
#' @param nrepeat The number of repetitions for cross-validation.For a thorough tuning step, the following code should be repeated 10 - 50 times and the error rate is averaged across the runs.
#' @param dist only applies to an object inheriting from "plsda" or "splsda" to evaluate the classification performance of the model. Should be a subset of "max.dist", "centroids.dist", "mahalanobis.dist".
#' @param measure Metrics for evaluating classification performance.
#' @param validation character. What kind of (internal) validation to use, matching one of "Mfold" or "loo". Default is "Mfold".
#' @param invisible If true, only the final model is output, otherwise the final model and classification verification results are output.
#' @param nCore multi-core computing.
#' @import mixOmics
#' @return splsda classification model or a list containing the classification model and 10X cross-classification validation results
#' @export

sPLSDA <- function(X,Y,ncomp,KeepX,folds = 10,nrepeat = 10,dist = 'max.dist',measure = "BER",validation = 'Mfold',invisible = T,nCore=1){
  # X is a matrix with rows representing cells/spots and columns representing genes, while Y is a character vector representing labels.
  if (!(is.matrix(X) || is.data.frame(X) || is(X, "sparseMatrix"))) {
    stop("X must a matrix(dcgmatrix) or data.frame")
  }
  # remove rare cells/spots
  non_rare_cell <- names(table(Y)[table(Y) > 10])
  X <- X[Y %in% non_rare_cell,]
  Y <- Y[Y %in% non_rare_cell]

  #Determine optimal component quantities and composition based on metrics and 10X cross-validation
  tune.splsda <-tune.splsda(X,Y,ncomp = ncomp,
                            validation = validation, folds = folds, nrepeat = nrepeat,
                            dist = dist,measure = measure,
                            test.keepX = KeepX,cpus = nCore)
  ncomp <- tune.splsda$choice.ncomp$ncomp
  select.keepX <- tune.splsda$choice.keepX
  splsda <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX)

  #Feature selection based on 10X cross-validation.
  perf.splsda<- perf(splsda, folds = folds, validation = validation,
                     dist = dist, progressBar = FALSE, nrepeat = nrepeat)
  if (!invisible) {
    splsda <- list(tune.splsda,splsda,perf.splsda)
  }
  return(splsda)
}


#' @title Screening for stable characteristic genes.
#'
#' @param perf Model cross-validation results output by sPLSDA.
#' @param stable Frequency of feature genes being selected in cross-validation.
#' @param org Corresponding species whole genome annotation.
#'
#' @return A matrix containing characteristic genes and corresponding stabilities.
#' @export
#' @import magrittr
#' @import clusterProfiler
featureSelect <- function(perf, stable = 0.9, org = NULL) {
  if (!inherits(perf, "perf")) {
    stop('Input must be the output result of sPLSDA (perf object)')
  }
  if (is.null(org)) {
    stop('Gene name conversion is required, please select the correct Genome annotation file!')
  }
  identify_gene_format <- function(gene_name) {
    if (grepl("^ENSG[0-9]+", gene_name)) {
      return("ENSYMBL")
    } else if (grepl("^[A-Za-z0-9]+$", gene_name)) {
      return("SYMBOL")
    } else if (grepl("^[0-9]+$", gene_name)) {
      return("ENTREZID")
    } else {
      stop("Unknown GeneID fotmat !")
    }
  }
  genetype <- identify_gene_format(names(perf[["features"]][["stable"]][[1]][1]))
  genelist = data.frame()
  for (i in names(perf$predict)) {
    if (genetype != 'ENTREZID') {
      Gene_ID <-
        bitr(names(perf[["features"]][["stable"]][[i]])[perf[["features"]][["stable"]][[i]] >= stable],
             fromType = genetype,
             toType = "ENTREZID",
             OrgDb = org)
    } else{
      Gene_ID <- names(perf[["features"]][["stable"]][[i]])
    }
    if (length(Gene_ID) == 0) {
      genelist = genelist
    } else{
      tmp = data.frame(
        Gene_ID,
        group = i,
        stable = as.numeric(perf[["features"]][["stable"]][[i]][Gene_ID[[genetype]]])
      )
      genelist = rbind(genelist, tmp)
    }
  }
  return(genelist)
}










### 可视化 ####

#' @title Plot the spatial distribution of component scores.
#'
#' @param sPLS A sPLS-DA model
#' @param comps A character vector containing the components to be involved in the plot, which must match the column names format of comp_score.
#' @param spatial A data frame, containing spatial coordinates (x, y), and labels (label).
#' @param size The size of the points.
#' @import ggplot2
#' @import cowplot
#' @import dplyr
#' @import RColorBrewer
#' @return A list of ggplot objects, and prints a combined plot.
#' @export
#'
#' @examples
#' set.seed(123)
#' sPLS <- list(variates = list(X = data.frame(comp1 = runif(10), comp2 = runif(10), comp3 = runif(10))))
#' spatial <- data.frame(x = runif(10), y = runif(10))
#' PlotComppos(sPLS, comps = 'ALL', spatial = spatial, size = 0.45)
#'
PlotComppos <- function(sPLS,comps = 'ALL',spatial,size=0.45){
  comp_score <- as.data.frame(sPLS$variates$X)
  if(comps == 'ALL') {comps <- colnames(comp_score)}
  else if(length(setdiff(comps,colnames(comp_score))) != 0){
    stop("Invalid index! Please check the component names and enter the correct component.")}
  plot_list <- list()  # Initialize an empty list to store plots
  for (i in comps) {
    p <- ggplot(data = comp_score) +  # Start a ggplot object with comp_score data
      geom_rect(aes(xmin = spatial$x - size,  # Define rectangles' x and y boundaries
                    xmax = spatial$x + size,
                    ymin = spatial$y - size,
                    ymax = spatial$y + size,
                    fill = !!sym(i))) +  # Fill rectangles with component data
      theme_bw(base_size = 14) +  # Set theme with base font size
      theme(panel.grid.major = element_blank(),  # Remove major grid lines
            panel.grid.minor = element_blank()) +  # Remove minor grid lines
      labs(title = i) +  # Set plot title to the component name
      theme(axis.ticks = element_blank(),  # Remove axis ticks
            axis.line = element_blank(),  # Remove axis lines
            axis.text = element_blank(),  # Remove axis text
            axis.title = element_blank(),  # Remove axis titles
            legend.text = element_text(size = 10),  # Set legend text size
            panel.border = element_rect(color = 'grey')) +  # Set panel border color
      theme(legend.position = "right", legend.title = element_blank()) +  # Position legend on the right
      guides(color = 'none') +  # Remove color guide
      theme(plot.title = element_text(hjust = 0.5)) +  # Center plot title
      scale_fill_viridis_c() +  # Apply Viridis color scale
      theme(text = element_text(family = "serif"))  # Set font family to serif
    plot_list[[i]] <- p  # Store the plot in the list with the component name as key
  }
  combined_plot <- plot_grid(plotlist = plot_list, align = "hv")
  return(combined_plot)
}






#' @title Visualization of clusters.
#'
#' @param Seurat_Obj Seurat object, the result of Seurat clustering, and the colored label is the corresponding annotation information.
#' @param subtype The annotation category to focus on, if not specified, select all.
#' @param type The label information used for coloring needs to be stored in the colData of the `seuratObject`.
#' @param size Scatter plot point size.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @return A ggplot object
#'
#'
PlotCsub <- function(Seurat_Obj,type = c('celltype','layer','seurat_clusters'),subtype = NULL,size=1){
  library(RColorBrewer)
  if (! type %in% colnames(Seurat_Obj@meta.data)) {
    stop('Cell or spot annotations are not included in the metadata of the Seurat object !')
  }
  if (!'umap' %in% names(Seurat_Obj@reductions)) {
    stop('Umap embeddings cannot be found, please perform dimensionality reduction first !')
  }
    df <- data.frame(Embeddings(object = Seurat_Obj,reduction = 'umap'))
    if (type == 'seurat_clusters') {
      df$label <- Seurat_Obj$seurat_clusters
    }else{df$label <- Seurat_Obj[[type]]}
    colnames(df) <- c('X','Y','label')
    if (is.null(subtype)) {
      subtype <- unique(df$label[[1]])
    }
    df$color <- "others"
    n_colors <- length(subtype)
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
    new_color <- getPalette(n_colors)
    color_mapping <- new_color[1:length(subtype)]
    names(color_mapping) <- subtype
    color_mapping <- c(color_mapping,"others"="grey")

    for (i in 1:length(subtype)) {
      df$color[df$label==subtype[i]] <- subtype[i]
    }
    p <- ggplot(df, aes(x=.data$X, y=.data$Y, color = .data$color)) +
      geom_point(size=size) +
      scale_color_manual(values = color_mapping) +
      guides(color = guide_legend(title = type,override.aes = list(size = 10))) +
      theme_bw() +
      theme(panel.grid.major=element_line(colour=NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank())
  return(p)
}


#' @title Label positioning after dimensionality reduction
#'
#' @param umap_df A dataframe containing UMAP coordinates.
#' @param Seurat_Obj A Seurat object containing clustering information or other metadata.
#' @param size Size of the points in the plot.
#' @param type The type of label to use for coloring points. Default is "label".
#' @param subtype Optional argument specifying which subtypes (clusters) to highlight.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @return A ggplot object

PlotCompumap <- function(umap_df,Seurat_Obj,size,type="label",subtype = NULL){
  library(RColorBrewer)
  if (type=='seurat_clusters') {
    umap_df$label <- Seurat_Obj$seurat_clusters
  }else{
    umap_df$label <- Seurat_Obj[[type]][[1]]
  }
  colnames(umap_df) <- c('X','Y',type)
  cluster <- unique(umap_df[[type]])

  if (is.null(subtype)) {
    subtype <- cluster
  }

  umap_df$color <- "others"
  n_colors <- length(subtype)
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  new_color <- getPalette(n_colors)
  color_mapping <- new_color[1:length(subtype)]
  names(color_mapping) <- subtype
  color_mapping <- c(color_mapping,"others"="grey")

  for (i in 1:length(subtype)) {
    umap_df$color[umap_df[[type]]==as.character(subtype[i])] <- as.character(subtype[i])
  }

  p <- ggplot(umap_df, aes(x=.data$X, y=.data$Y, color = .data$color)) +
    geom_point(size=size) +
    scale_color_manual(values = color_mapping) +
    guides(color = guide_legend(title = type,override.aes = list(size = 10))) +
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())
  return(p)
}


#' @title Contribution of principal components
#'
#' @param GO_res Results of GO or other enrichment analysis.
#' @param sPLS The final discriminant model of `sPLSDA`.
#' @param gene_list Gene set of interest, if not specified, selects all genes.
#' @param nterm The number of terms that need to be displayed for enrichment analysis.
#' @param OrgDb Corresponding species whole genome annotation
#' @param comp Principal components to be displayed
#'
#' @import ggplot2
#' @import clusterProfiler
#' @import dplyr
#' @return A bar chart
#' @export
#'

PlotLoad <- function(GO_res,sPLS,gene_list=NULL,nterm=NULL,comp,org='org.Hs.eg.db'){

  if("rfe" %in% class(sPLS)) {
    model_pls <- sPLS[["fit"]][["finalModel"]]
    temp_loadings <- loadings(model_pls)
    colname <- NULL
    loadings_df <- NULL
    for (i in 1:ncol(temp_loadings)) {
      colname <- c(colname,paste("comp",i,sep = ""))
      loadings_df <- cbind(loadings_df,temp_loadings[,i])
    }
    loadings_df <- data.frame(loadings_df)
    colnames(loadings_df) <- colname
    feature <- sPLS$optVariables
  }else{
    loadings_df <- data.frame(sPLS[["loadings"]][["X"]])
    feature <- gene_list$ENTREZID[gene_list$group == comp]
  }

  df <- data.frame(gene = feature,
                   loadings = loadings_df[feature,comp])
  sorted_df <- df[order(df$loadings), ]

  # 提取gene所在的term信息
  if (is.null(nterm)) {
    GO_res <- data.frame(GO_res)
  }else{
    GO_res <- na.omit(data.frame(GO_res[1:nterm]))
  }


  if ('core_enrichment' %in% colnames(GO_res)) {
    GO_set <- str_split(string = GO_res$core_enrichment,pattern = '/',simplify = F)
  }else{
    GO_set <- str_split(string = GO_res$geneID,pattern = '/',simplify = F)
  }


  names(GO_set) <- GO_res$Description
  term_to <- do.call(rbind,
                     lapply(names(GO_set),
                            function(name) {data.frame(gene = GO_set[[name]], term = name,stringsAsFactors = T)})
  )


  if (!'ONTOLOGY' %in% colnames(GO_res)) {
    symbol <- bitr(term_to$gene, fromType="ENTREZID", toType="SYMBOL", OrgDb=org)
    names(symbol) <- c( "gene","SYMBOL")
    term_to <- left_join(term_to,symbol,"gene")
    names(term_to) <- c("ENTREZID","term","gene")
  }

  # combind
  gene_loadings_term <- left_join(x = sorted_df,y = term_to,by = "gene")

  # barplot
  cluster <- unique(gene_loadings_term$term)
  n_colors <- length(cluster)
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  new_color <- getPalette(n_colors)
  color_mapping <- new_color[1:length(cluster)]
  names(color_mapping) <- cluster


  p <- ggplot(na.omit(gene_loadings_term), aes(x = loadings, y = reorder(gene,loadings), fill = term)) +
    geom_bar(stat = "identity", position = "identity") +
    labs(title = comp, x = "loadings", y = "gene") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 7)) +
    scale_fill_manual(values = color_mapping) +
    theme_minimal()
  return(p)
}
