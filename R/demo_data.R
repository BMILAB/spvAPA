#' Gene Expression Matrix for ST11GEM
#'
#' This dataset contains a gene expression matrix where rows represent genes
#' and columns represent cells for the ST11GEM experiment.
#'
#' @format A numeric matrix:
#' \describe{
#'   \item{Rows}{Genes, where each row corresponds to a specific gene.}
#'   \item{Columns}{Spots, where each column corresponds to a specific spot.}
#' }
#' @source Provided by the user for the ST11GEM experiment.
#' @examples
#' data(ST11GEM)
#' dim(ST11GEM)  # View dimensions of the matrix
#' head(ST11GEM[, 1:5])  # View first 5 cells for the first few genes
"ST11GEM"

#' Gene Expression Matrix for ST11RUD
#'
#' This dataset contains a RUD matrix where rows represent genes
#' and columns represent cells for the ST11RUDa experiment.
#'
#' @format A numeric matrix:
#' \describe{
#'   \item{Rows}{Genes, where each row corresponds to a specific gene.}
#'   \item{Columns}{Cells, where each column corresponds to a specific spot.}
#' }
#' @source Provided by the user for the ST11RUD experiment.
#' @examples
#' data(ST11RUD)
#' dim(ST11RUD)  # View dimensions of the matrix
#' head(ST11RUD[, 1:5])  # View first 5 cells for the first few genes
"ST11RUD"

#' Metadata for ST11label
#'
#' This dataset contains metadata for the ST11 experiments, including
#' coordinates and annotation labels for cells.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{x}{Numeric. The x-coordinate of the cell.}
#'   \item{y}{Numeric. The y-coordinate of the cell.}
#'   \item{label}{Character. The annotation label for the cell (e.g., cell type).}
#' }
#' @source Provided by the user for the ST11 experiments.
#' @examples
#' data(ST11label)
#' head(ST11label)  # View the first few rows of metadata
#' table(ST11label$label)  # Summarize cell type annotations
"ST11label"
