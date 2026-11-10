# Namespace declarations for base graphics/statistics helpers. Data columns in
# dplyr verbs and ggplot2 aesthetics are referenced through the rlang `.data`
# pronoun, so no global variables are declared.

#' @importFrom data.table :=
#' @importFrom grDevices as.raster chull dev.off pdf
#' @importFrom graphics abline arrows par points text
#' @importFrom methods new
#' @importFrom rlang .data
#' @importFrom stats as.dendrogram as.dist cor density dnorm hclust median
#'   model.matrix na.omit p.adjust prcomp quantile setNames
#' @importFrom utils tail
#' @noRd
NULL
