############ Adapted functions from Cardinal for binning matrix objects ############

#' Spectral binning of intensity values stored in a Matrix object, converted
#' from matter.
#'
#' @param matrix Numeric matrix of features by pixels.
#' @param ref Positive reference mass-to-charge (m/z) values for the output bins.
#' @param index Numeric m/z values, one per input matrix row.
#' @param method A character vector specifying the binning method. Options
#'   include `"sum"`, `"mean"`, `"max"`, `"min"`. If not specified default
#'   method used is "sum".
#' @param tolerance Positive mass tolerance in Da. A value named `relative`
#'   is interpreted as a fraction of each reference mass (for example,
#'   `c(relative = 5e-6)` for 5 ppm).
#'
#' @return Matrix object containing the binned intensity values matching the
#'   provided reference list.
#' @export
#'
#' @examples
#' intensity <- matrix(1:12, nrow = 3)
#' spectralBinning(intensity,
#'     ref = c(100, 200), index = c(100, 100.01, 200),
#'     method = "sum", tolerance = c(absolute = 0.02)
#' )
spectralBinning <- function(
    matrix,
    ref,
    index,
    method = c("sum", "mean", "max", "min"),
    tolerance
) {
    method <- match.arg(method)
    if (!is.matrix(matrix) || !is.numeric(matrix)) {
        stop("'matrix' must be a numeric matrix.", call. = FALSE)
    }
    if (!is.numeric(ref) || !length(ref) || any(!is.finite(ref) | ref <= 0)) {
        stop("'ref' must contain positive finite reference masses.", call. = FALSE)
    }
    if (!is.numeric(index) || length(index) != nrow(matrix) ||
        any(!is.finite(index) | index <= 0)) {
        stop("'index' must contain one positive finite mass per matrix row.", call. = FALSE)
    }
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance <= 0 ||
        (!is.null(names(tolerance)) && !names(tolerance) %in% c("absolute", "relative"))) {
        stop("'tolerance' must be one positive finite number, optionally named 'absolute' or 'relative'.", call. = FALSE)
    }
    widths <- if (identical(names(tolerance), "relative")) ref * tolerance else rep(tolerance, length(ref))
    binned_matrix <- matrix(0, nrow = length(ref), ncol = ncol(matrix),
        dimnames = list(as.character(ref), colnames(matrix)))
    for (i in seq_along(ref)) {
        in_bin <- which(abs(index - ref[i]) <= widths[i])
        if (length(in_bin)) {
            values <- matrix[in_bin, , drop = FALSE]
            binned_matrix[i, ] <- switch(method,
                sum = colSums(values), mean = colMeans(values),
                max = apply(values, 2, max), min = apply(values, 2, min))
        }
    }
    binned_matrix
}


#' @noRd
.binSpaMTPMatrix <- function(
    data,
    resolution,
    units,
    assay,
    assayName,
    method
) {
    units <- match.arg(units, c("ppm", "mz"))
    method <- match.arg(method, c("sum", "mean", "max", "min"))
    if (!is.numeric(resolution) || length(resolution) != 1L ||
        !is.finite(resolution) || resolution <= 0) {
        stop("resolution must be one positive finite number.", call. = FALSE)
    }
    expression <- .assayData(data, assay, assayName)
    mz <- .massValues(.featureMetadata(data, assay), rownames(expression))
    if (!nrow(expression) || length(mz) != nrow(expression)) {
        stop("Binning needs an m/z value for every feature.", call. = FALSE)
    }
    origin <- min(mz)
    if (units == "ppm" && origin <= 0) {
        stop(
            "Relative ppm binning requires positive m/z values.",
            call. = FALSE
        )
    }
    # Assign each feature exactly once; exact midpoints go to the upper bin.
    # Work on the log-mass axis for constant relative (ppm) spacing.
    position <- if (units == "ppm") {
        log(mz / origin) / log1p(resolution * 1e-6)
    } else {
        (mz - origin) / resolution
    }
    if (any(!is.finite(position)) || any(position > 2^52)) {
        stop(
            "resolution is too small for numerically stable binning.",
            call. = FALSE
        )
    }
    bin <- floor(position + 0.5)
    occupied <- sort(unique(bin))
    group <- match(bin, occupied)
    reference <- if (units == "ppm") {
        origin * exp(occupied * log1p(resolution * 1e-6))
    } else {
        origin + occupied * resolution
    }
    membership <- Matrix::sparseMatrix(
        i = group, j = seq_along(group), x = 1,
        dims = c(length(occupied), nrow(expression))
    )
    if (method %in% c("sum", "mean")) {
        result <- membership %*% expression
        if (method == "mean") {
            result <- Matrix::Diagonal(x = 1 / tabulate(group, length(occupied))) %*% result
        }
    } else {
        result <- t(vapply(seq_along(occupied), function(index) {
            apply(
                as.matrix(expression[group == index, , drop = FALSE]), 2L,
                if (method == "max") max else min
            )
        }, numeric(ncol(expression))))
    }
    dimnames(result) <- list(as.character(reference), colnames(expression))
    result
}

#' Bin aligned spatial metabolomics features
#'
#' `MSImagingArrays` or `MSImagingExperiment` input delegates binning to
#' Cardinal and is converted
#' to SpatialExperiment only after a common feature axis exists.
#' SpatialExperiment input returns a new binned SpatialExperiment.
#' For aligned matrices, bins start at the minimum observed m/z and use
#' linear (`mz`) or log-mass (`ppm`) spacing. Each feature is assigned once
#' to its nearest bin on that axis (midpoints go up); empty bins are omitted.
#' Thus sum binning conserves each pixel's total intensity. Cardinal input
#' uses Cardinal's reference-window binning convention instead.
#'
#' @param data A `SpatialExperiment`, Cardinal `MSImagingArrays` or
#'   `MSImagingExperiment`.
#' @param resolution Positive bin resolution.
#' @param units Resolution units, either `"ppm"` or `"mz"`.
#' @param assay Primary MSI experiment name; alternative modalities cannot be
#'   binned.
#' @param assayName Matrix layer or SpatialExperiment assay to bin.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param method Aggregation method.
#' @param returnOnlyMtx Return only the binned matrix.
#'
#' @return A binned `SpatialExperiment` or matrix.
#' @export
#' @examples
#' methods::showMethods("binSpaMTP")
methods::setGeneric(
    "binSpaMTP",
    function(
        data,
                     resolution,
                     units = "ppm",
                     assay = "main",
                     assayName = "counts",
                     method = c("sum", "mean", "max", "min"),
                     returnOnlyMtx = FALSE
    ) {
        methods::standardGeneric("binSpaMTP")
    }
)

#' @rdname binSpaMTP
#' @export
methods::setMethod(
    "binSpaMTP",
    "SpatialExperiment",
    function(
        data,
                     resolution,
                     units = "ppm",
                     assay = "main",
                     assayName = "counts",
                     method = c("sum", "mean", "max", "min"),
                     returnOnlyMtx = FALSE
    ) {
        method <- match.arg(method)
        if (assay %in% SingleCellExperiment::altExpNames(data)) {
            stop(
                "Bin the primary MSI experiment, not an alternative modality.",
                call. = FALSE
            )
        }
        matrix <- .binSpaMTPMatrix(
            data,
            resolution,
            units,
            assay,
            assayName,
            method
        )
        if (isTRUE(returnOnlyMtx)) {
            return(matrix)
        }
        mz <- suppressWarnings(as.numeric(sub("^mz-", "", rownames(matrix))))
        if (any(!is.finite(mz))) {
            mz <- .massValues(.featureMetadata(data, assay), rownames(matrix))
        }
        featureNames <- make.unique(
            paste0("mz-", format(mz, digits = 15, trim = TRUE))
        )
        rownames(matrix) <- featureNames
        object <- SpatialExperiment::SpatialExperiment(
            assays = list(counts = matrix),
            rowData = S4Vectors::DataFrame(
                mz = mz,
                raw_mz = mz,
                mz_names = featureNames,
                row.names = featureNames
            ),
            colData = SummarizedExperiment::colData(data),
            spatialCoords = SpatialExperiment::spatialCoords(data),
            imgData = SpatialExperiment::imgData(data),
            metadata = S4Vectors::metadata(data)
        )
        for (name in SingleCellExperiment::altExpNames(data)) {
            SingleCellExperiment::altExp(object, name) <-
                SingleCellExperiment::altExp(data, name)
        }
        S4Vectors::metadata(object)$spamtp_binning <- list(
            resolution = resolution,
            units = units,
            method = method
        )
        # Peak-level results and graphs refer to the pre-binning representation.
        for (name in c(
            "mz_annotation", "db_3", "moransi", "spatialGraphs",
            "reductionLoadings", "spamtp_normalization", "spamtp_integration"
        )) {
            S4Vectors::metadata(object)[[name]] <- NULL
        }
        object
    }
)

#' @rdname binSpaMTP
#' @export
methods::setMethod(
    "binSpaMTP",
    "MSImagingExperiment",
    function(
        data,
                     resolution,
                     units = "ppm",
                     assay = "main",
                     assayName = "counts",
                     method = c("sum", "mean", "max", "min"),
                     returnOnlyMtx = FALSE
    ) {
        .binCardinal(
            data,
            resolution,
            units,
            match.arg(method),
            returnOnlyMtx
        )
    }
)

.binCardinal <- function(data, resolution, units, method, returnOnlyMtx) {
    binned <- .binCardinalData(data, resolution, units, method)
    if (isTRUE(returnOnlyMtx)) {
        return(Cardinal::spectra(binned))
    }
    cardinalToSpatialExperiment(binned)
}

.binCardinalData <- function(
    data, resolution, units, method,
    massRange = NULL, verbose = NULL
) {
    if (length(resolution) != 1L || !is.numeric(resolution) ||
        !is.finite(resolution) || resolution <= 0) {
        stop("resolution must be one positive finite number.", call. = FALSE)
    }
    data <- Cardinal::process(data, verbose = FALSE)
    arguments <- list(
        x = data, resolution = resolution,
        units = match.arg(units, c("ppm", "mz")), method = method
    )
    if (!is.null(massRange)) arguments$mass.range <- massRange
    if (!is.null(verbose)) arguments$verbose <- verbose
    if (methods::is(data, "MSImagingArrays") && is.null(massRange)) {
        # Cardinal's default reference is inferred from the first spectrum.
        # Inspect one mass vector at a time to cover all pixels without creating
        # a dense intensity matrix or discarding peaks unique to later spectra.
        ranges <- vapply(Cardinal::mz(data), function(values) {
            if (!length(values) || any(!is.finite(values))) {
                stop(
                    "Every raw spectrum needs finite m/z values before binning.",
                    call. = FALSE
                )
            }
            range(values)
        }, numeric(2))
        arguments$mass.range <- c(floor(min(ranges)), ceiling(max(ranges)))
    }
    do.call(Cardinal::bin, arguments)
}

#' @rdname binSpaMTP
#' @export
methods::setMethod(
    "binSpaMTP", "MSImagingArrays",
    function(
        data, resolution, units = "ppm", assay = "main", assayName = "counts",
                     method = c("sum", "mean", "max", "min"), returnOnlyMtx = FALSE
    ) {
        .binCardinal(
            data,
            resolution,
            units,
            match.arg(method),
            returnOnlyMtx
        )
    }
)

#' @rdname binSpaMTP
#' @export
methods::setMethod(
    "binSpaMTP",
    "ANY",
    function(
        data,
                     resolution,
                     units = "ppm",
                     assay = "main",
                     assayName = "counts",
                     method = c("sum", "mean", "max", "min"),
                     returnOnlyMtx = FALSE
    ) {
        .requireExperiment(data, "SpatialExperiment")
    }
)
