#' Pathway-network reaction styles
#'
#' A small lookup table mapping SpaMTP interaction style codes to visual
#' properties used
#' by the interactive pathway-network viewer.
#'
#' @format A data frame with 18 rows and 5 variables:
#' \describe{
#'   \item{reaction_type}{Stable SpaMTP style code; not a database-specific factor index.}
#'   \item{reaction_name}{Human-readable reaction type.}
#'   \item{linetype}{Line style.}
#'   \item{arrowhead}{Arrowhead style.}
#'   \item{colour}{Edge colour.}
#' }
#' @return A data frame of pathway-network reaction styles.
#' @source SpaMTP-maintained display vocabulary. The source interaction labels
#'   are mapped to these stable style codes in
#'   `data-raw/pathway_interaction_utils.R` in the SpaMTPdb source repository.
#'   This table defines visual styles rather than biological measurements.
#' @keywords datasets
"reaction_type"
