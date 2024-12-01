## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4,
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
library(spvAPA)

## ----eval=TRUE----------------------------------------------------------------
data("ST11GEM",package = "spvAPA")
data("ST11RUD",package = "spvAPA")
data("ST11label",package = "spvAPA")
dim(ST11GEM)
dim(ST11RUD)

## -----------------------------------------------------------------------------
RUD_imp <- WNNImpute(gene = ST11GEM,APA = ST11RUD,k = 15,init = F,is.weight = F)

## -----------------------------------------------------------------------------
RUD_imp_t <- t(RUD_imp)
layer <- ST11label$label
list.keepX <- c(1:10,  seq(20, 150, 10))
RUD_spls <- sPLSDA(X = RUD_imp_t,Y = layer,ncomp = 4,KeepX = list.keepX,folds = 10,invisible = F,nCore = 12)

## -----------------------------------------------------------------------------
GEM_t <- t(ST11GEM)
layer <- ST11label$label
list.keepX <- c(1:10,  seq(20, 150, 10))
GEM_spls <- sPLSDA(X = GEM_t,Y = layer,ncomp = 4,KeepX = list.keepX,folds = 10,invisible = F,nCore = 12)

## -----------------------------------------------------------------------------
library(org.Mm.eg.db)
RUD_perf <- RUD_spls[[3]]
featureSelect(perf = RUD_perf,stable = 0.9)

## -----------------------------------------------------------------------------
# Obtain the score matrix of components.
RUD_mod <- RUD_spls[[2]]
comp_score <- as.data.frame(RUD_mod$variates$X)

# Plot

library(cowplot)

comp_score_plot <- PlotComppos(sPLS = RUD_mod,label = ST11label,comps = "ALL",size = 1.5)
plot_grid(plotlist = comp_score_plot, align = "hv")

## -----------------------------------------------------------------------------
PlotCompumap(spls = RUD_mod,size = 1.7)

## -----------------------------------------------------------------------------
GEM_mod <- GEM_spls[[2]]
PlotCompumap(spls = RUD_mod,size = 1.7)

## -----------------------------------------------------------------------------
PlotCompumap(spls = list(RUD_mod,GEM_mod),size = 1.7)

