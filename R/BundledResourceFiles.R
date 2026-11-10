# Development assistance: see inst/COPYRIGHTS.

#' Provenance of bundled pathway corrections and spatial calibration
#'
#' Small support files are installed in `extdata`. They can be located with
#' `system.file("extdata", filename, package = "SpaMTP")`.
#'
#' Pathway interaction corrections and their provenance are distributed by
#' SpaMTPdb; see [SpaMTPdb::spaMTPdbResource()].
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
