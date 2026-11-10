# Development assistance: see inst/COPYRIGHTS.

#' Small offline inputs for the function examples
#'
#' Deterministic teaching data with eight mass features and eight paired genes
#' at 24 spatial pixels. Intensities, regions, gene identifiers and pathways
#' are synthetic. The three chemical formulae illustrate mass matching only;
#' the pathway memberships have no biological interpretation.
#' Preparation code is in `data-raw/create_spamtp_example.R` in the source
#' repository.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{experiment}{An annotated SpatialExperiment with counts, logcounts,
#'     spatial coordinates, region and replicate labels, a PCA embedding and
#'     paired gene counts in the transcriptome alternative experiment.}
#'   \item{metabolites}{Four chemical records with masses and SMILES features,
#'     including two isomers to illustrate ambiguous annotation.}
#'   \item{database}{A miniature RaMP-shaped resource bundle with three
#'     fictional pathways and explicitly supplied gene mappings.}
#'   \item{gene_reference}{Eight fictional HGNC-shaped identifier records.}
#' }
#' @return A list of example inputs.
#' @examples
#' data(spamtp_example)
#' plotSpatialFeature(spamtp_example$experiment, "ion1")
#' @keywords datasets
"spamtp_example"
