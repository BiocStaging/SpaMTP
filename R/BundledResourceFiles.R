# Development assistance: see inst/COPYRIGHTS.

#' Provenance of bundled pathway corrections and spatial calibration
#'
#' Small support files are installed in `extdata`. They can be located with
#' `system.file("extdata", filename, package = "SpaMTP")`.
#'
#' @section Pathway correction maps:
#' `ramp_kegg_interactions_v1.rds`, `ramp_reactome_interactions_v1.rds`,
#' `ramp_wikipathway_interactions_v1.rds` and `ramp_hmdb_interactions_v1.rds`
#' are SpaMTP-generated corrections for the archived RaMP 3.0.7 collections.
#' The historical `ramp_hmdb` name refers to SMPDB topology.
#' Each file stores an input SHA-256 fingerprint, source provenance and,
#' for each pathway edge table, row indices, direction, SpaMTP style codes
#' and source interaction labels. Repeated row indices restore parallel edges.
#' These maps contain no node identifiers or pathway graphs. A correction
#' applies only when the input collection matches its recorded fingerprint.
#'
#' The recipe is `data-raw/rebuild_pathway_interactions.R`, with
#' `data-raw/pathway_interaction_utils.R` and `data-raw/ramp_graph_utils.R`,
#' in the source repository. It uses graphite archive 19 (2022-11-01) and
#' pinned SpaMTPdb resources. `pathway_interactions_provenance.csv` records
#' each input filename, source URL and MD5 checksum;
#' `pathway_interactions_audit.csv` records the numbers of pathways, edge
#' tables and edges before and after correction for each collection.
#' Source graphs retain the terms of their upstream databases; correction
#' does not change those terms. See the installed `COPYRIGHTS` file.
#'
#' @section Mouse-brain spatial calibration:
#' `mouse_brain_scalefactors.json` contains four Space Ranger factors for
#' the V11L12-109_B1 sample distributed by Vicari et al. It preserves the
#' high- and low-resolution image scales, fiducial diameter and spot diameter
#' in full-resolution image coordinates. The spot diameter is 376.942 pixels.
#' These values calibrate only that archived sample; other samples need their
#' own measured scales. The source dataset is licensed CC BY-NC 3.0.
#'
#' @source Graphite pathway archive,
#'   \url{https://graphiteweb.bio.unipd.it/pathways/19/};
#'   RaMP snapshot used by the correction recipe,
#'   \doi{10.5281/zenodo.22045311};
#'   Vicari et al., Spatial Multimodal Analysis of Transcriptomes and
#'   Metabolomes in Tissues, \doi{10.17632/w7nw4km7xd.1}.
#' @seealso [loadSpaMTPDatabase()], [vignetteData], [reaction_type]
#' @return This documentation topic describes installed support files;
#'   it does not define a function.
#' @name bundledResourceFiles
NULL
