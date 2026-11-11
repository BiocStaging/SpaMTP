#### SpaMTP Cardinal SSC clustering segmentation ###############################################################################################################################################################################


#' Adds Cardinal ssc segmentation annotation to m/z count data object
#'
#' Adds spatial shrunken centroid (SSC) labels to
#' `Cardinal::pixelData(data)[["ssc"]]` through its public replacement accessor.
#' This function accepts Cardinal containers only; it does not modify S4 slots
#' directly or require Seurat.
#'
#' @param data A Cardinal Object containing the raw/binned m/z count data.
#' @param dataSsc A Cardinal Object containing the ssc segmentation results.
#'   Note: Cardinal's spatialShrunkenCentroids() must be run to generate this
#'   object.
#' @param resolution An integer (< v3.6.0) or string (>= v3.6.0) defining the
#'   ssc segmentation resolution to add to the raw Cardinal object.
#'
#' @returns A Cardinal Object containing the raw m/z counts and the relative ssc
#'   segmentation annotations per pixel.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' msi <- spatialExperimentToCardinal(spe)
#' segmentation <- withr::with_seed(1, Cardinal::spatialShrunkenCentroids(
#'     msi,
#'     r = 1, k = 2, s = 0, verbose = FALSE
#' ))
#' labelled <- addSSCAnnotation(msi, segmentation, resolution = 1)
#' table(Cardinal::pixelData(labelled)$ssc)
addSSCAnnotation <- function(data, dataSsc, resolution) {
    data_bin <- data

    if (check_cardinal_version()) {
        if (!methods::is(dataSsc, "ResultsList")) {
            Cardinal::pixelData(data_bin)[["ssc"]] <- dataSsc$class
        } else {
            message("Getting cluster segments for resolution ", resolution)

            Cardinal::pixelData(data_bin)[["ssc"]] <- dataSsc[[resolution]]$class
        }
    } else {
        message("Getting cluster segments for resolution (s) = ", resolution)
        cluster_idx <- which(Cardinal::modelData(dataSsc)[["s"]] == resolution)

        classes <- Cardinal::resultData(dataSsc)[[cluster_idx]][[1]]
        pixel_data <- Cardinal::pixelData(data_bin)


        pixel_data[["ssc"]] <- classes

        Cardinal::pixelData(data_bin) <- pixel_data
    }

    return(data_bin)
}

########################################################################################################################################################################################################################
