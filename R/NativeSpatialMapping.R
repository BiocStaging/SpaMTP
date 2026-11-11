.pixelPolygons <- function(coordinates, width, crs = NA_character_) {
    if (is.null(width)) {
        widths <- lapply(
            split(seq_len(nrow(coordinates)), coordinates$sample_id),
            function(index) {
                if (length(index) < 2L) {
                    return(numeric())
                }
                FNN::get.knn(
                    as.matrix(coordinates[index, c("x", "y")]),
                    k = 1
                )$nn.dist[, 1]
            }
        )
        distances <- unlist(widths)
        width <- stats::median(distances[distances > 0])
    }
    if (length(width) != 1L || !is.finite(width) || width <= 0) {
        stop(
            "Supply a positive smPixelWidth in spatial coordinate units.",
            call. = FALSE
        )
    }
    geometry <- lapply(seq_len(nrow(coordinates)), function(index) {
        square <- getSquareCoordinates(
            coordinates$x[index], coordinates$y[index],
            width, coordinates$cell[index]
        )
        sf::st_polygon(list(as.matrix(square[, c("X", "Y")])))
    })
    sf::st_sf(coordinates, geometry = sf::st_sfc(geometry, crs = crs))
}

.mappingWeights <- function(
    source, target, width, highResolution,
    radius, polygons, threshold
) {
    if (!length(intersect(source$sample_id, target$sample_id))) {
        stop("Source and target sample_id values do not match.", call. = FALSE)
    }
    crs <- if (is.null(polygons)) NA_character_ else sf::st_crs(polygons)
    sourcePolygons <- .pixelPolygons(source, width, crs)
    if (highResolution) {
        if (!is.null(polygons)) {
            stop("stPolygons is used for area-overlap mapping; set stHires = FALSE.",
                call. = FALSE
            )
        }
        targetGeometry <- sf::st_as_sf(target, coords = c("x", "y"), crs = crs)
    } else if (!is.null(polygons)) {
        if (!inherits(polygons, "sf") || !"cell" %in% colnames(polygons) ||
            anyDuplicated(polygons$cell) || !setequal(polygons$cell, target$cell)) {
            stop("stPolygons must be sf polygons with cell IDs matching ST.data.",
                call. = FALSE
            )
        }
        if (isTRUE(sf::st_is_longlat(polygons))) {
            stop("Use planar spatial coordinates and polygons, not longitude/latitude.",
                call. = FALSE
            )
        }
        targetGeometry <- polygons[match(target$cell, polygons$cell), ]
        if (any(!sf::st_is_valid(targetGeometry)) || any(sf::st_is_empty(targetGeometry))) {
            stop(
                "stPolygons must contain valid, non-empty polygons.",
                call. = FALSE
            )
        }
        if (any(!as.character(sf::st_geometry_type(targetGeometry)) %in%
            c("POLYGON", "MULTIPOLYGON"))) {
            stop("stPolygons must contain polygon geometries.", call. = FALSE)
        }
    } else {
        if (is.null(radius) || !length(radius) %in% c(1L, nrow(target)) ||
            any(!is.finite(radius)) || any(radius <= 0)) {
            stop(
                "Supply stRadius or stPolygons for area-overlap mapping.",
                call. = FALSE
            )
        }
        targetGeometry <- sf::st_buffer(
            sf::st_as_sf(target, coords = c("x", "y"), crs = crs),
            dist = rep(radius, length.out = nrow(target))
        )
    }
    candidates <- sf::st_intersects(targetGeometry, sourcePolygons)
    matches <- lapply(seq_along(candidates), function(index) {
        selected <- candidates[[index]]
        selected <- selected[source$sample_id[selected] == target$sample_id[index]]
        if (!length(selected)) {
            return(integer())
        }
        if (highResolution) {
            distance <- (source$x[selected] - target$x[index])^2 +
                (source$y[selected] - target$y[index])^2
            return(selected[order(distance, source$cell[selected])][1L])
        }
        targetArea <- as.numeric(sf::st_area(targetGeometry[index, ]))
        if (!is.finite(targetArea) || targetArea <= 0) {
            stop("Target polygons must have positive areas.", call. = FALSE)
        }
        coverage <- vapply(selected, function(pixel) {
            intersection <- suppressWarnings(sf::st_intersection(
                sf::st_geometry(targetGeometry[index, ]),
                sf::st_geometry(sourcePolygons[pixel, ])
            ))
            sum(as.numeric(sf::st_area(intersection))) / targetArea
        }, numeric(1))
        selected[coverage > 0 & coverage >= threshold]
    })
    count <- lengths(matches)
    weights <- Matrix::sparseMatrix(
        i = as.integer(unlist(matches, use.names = FALSE)),
        j = rep(seq_along(matches), count),
        x = rep(1 / pmax(count, 1L), count),
        dims = c(nrow(source), nrow(target)), dimnames = list(source$cell, target$cell)
    )
    list(weights = weights, matches = matches)
}

#' Map metabolite intensities to transcriptomic pixels or cells
#'
#' Both inputs must already use the same planar coordinate units and matching
#' sample_id values. Mapping never links different samples. `stHires` uses target
#' centroids inside square MSI pixels; ties use the nearest MSI centre then
#' pixel ID. Otherwise each MSI pixel must cover `overlapThreshold` of the target
#' spot/polygon area. Matching MSI intensities are averaged with equal weights.
#'
#' @param smData,stData SpatialExperiment objects with aligned coordinates.
#' @param stHires Use point-in-pixel mapping for high-resolution target cells.
#' @param smAssay,stAssay Primary (`"main"`) or alternative experiment names
#'   in the source and target containers; see [experimentAccess]. These select
#'   modalities, not individual expression matrices.
#' @param smPixelWidth Square MSI pixel width in coordinate units. NULL
#'   estimates the median positive nearest-neighbour distance within samples.
#' @param overlapThreshold Minimum fraction of target area covered by a pixel.
#' @param annotations Preserve MSI rowData and annotation provenance.
#' @param addMetadata Include summaries of source pixel metadata in colData.
#' @param mergeUniqueMetadata Remove duplicate values in source summaries.
#' @param mapData Map all source assays; FALSE maps counts only.
#' @param newSptAssay Name of the transcriptome alternative experiment.
#' @param newSpmAssay Name of the primary metabolomics experiment.
#' @param verbose Show progress messages.
#' @param stRadius Target spot radius, scalar or one per target pixel, in
#'   spatial coordinate units; alternatively supply ST.polygons.
#' @param stPolygons An sf polygon object with a cell column matching ST.data.
#' @param dropUnmapped Remove targets without matching MSI pixels. FALSE keeps
#'   zero-filled intensities and marks these pixels with mapped = FALSE.
#' @return A SpatialExperiment with MSI assays on the target coordinates, target
#' transcriptomes in altExp, target images and metadata, and mapping provenance.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' mapped <- mapSpatialOmics(spe, spe,
#'     stHires = TRUE, smPixelWidth = 1,
#'     stAssay = "transcriptome"
#' )
#' table(mapped$mapped)
mapSpatialOmics <- function(
    smData, stData, stHires = FALSE, smAssay = "main", stAssay = "main",
    smPixelWidth = NULL, overlapThreshold = 0.2, annotations = TRUE,
    addMetadata = TRUE, mergeUniqueMetadata = TRUE, mapData = FALSE,
    newSptAssay = "transcriptome", newSpmAssay = "main", verbose = FALSE,
    stRadius = NULL, stPolygons = NULL, dropUnmapped = FALSE
) {
    smData <- .nativeSpatialObject(smData)
    stData <- .nativeSpatialObject(stData)
    source <- .nativeCoordinates(smData)
    target <- .nativeCoordinates(stData)
    if (length(overlapThreshold) != 1L || !is.finite(overlapThreshold) ||
        overlapThreshold < 0 || overlapThreshold > 1) {
        stop("overlapThreshold must be between zero and one.", call. = FALSE)
    }
    mapping <- .mappingWeights(
        source, target, smPixelWidth, stHires,
        stRadius, stPolygons, overlapThreshold
    )
    experiment <- .experimentForAssay(smData, smAssay)
    selected <- if (mapData) SummarizedExperiment::assayNames(experiment) else "counts"
    if (!length(selected) || !all(selected %in% SummarizedExperiment::assayNames(experiment))) {
        stop(
            "Source assays are missing; provide counts or set mapData = TRUE.",
            call. = FALSE
        )
    }
    assays <- stats::setNames(lapply(selected, function(name) {
        SummarizedExperiment::assay(experiment, name) %*% mapping$weights
    }), selected)
    featureData <- SummarizedExperiment::rowData(experiment)
    if (!annotations) {
        featureData <- S4Vectors::DataFrame(
            mz = .massValues(featureData, rownames(experiment)),
            row.names = rownames(experiment)
        )
    }
    metadata <- SummarizedExperiment::colData(stData)
    metadata$SPM_pixels <- S4Vectors::SimpleList(lapply(
        mapping$matches,
        function(index) source$cell[index]
    ))
    metadata$n_source_pixels <- lengths(mapping$matches)
    metadata$mapped <- lengths(mapping$matches) > 0L
    if (addMetadata) {
        sourceMetadata <- SummarizedExperiment::colData(smData)
        for (name in colnames(sourceMetadata)) {
            column <- tail(
                make.unique(c(colnames(metadata), paste0("SM_", name))),
                1L
            )
            metadata[[column]] <- vapply(mapping$matches, function(index) {
                values <- as.character(sourceMetadata[[name]][index])
                if (mergeUniqueMetadata) values <- unique(values)
                paste(values, collapse = ", ")
            }, character(1))
        }
    }
    output <- SpatialExperiment::SpatialExperiment(
        assays = assays, rowData = featureData, colData = metadata,
        spatialCoords = SpatialExperiment::spatialCoords(stData),
        imgData = SpatialExperiment::imgData(stData),
        metadata = S4Vectors::metadata(smData)
    )
    SingleCellExperiment::mainExpName(output) <- newSpmAssay
    SingleCellExperiment::altExp(output, newSptAssay) <-
        .experimentForAssay(stData, stAssay)
    S4Vectors::metadata(output)$spamtp_mapping <- list(
        method = if (stHires) "centroid within MSI pixel" else "target-area overlap",
        aggregation = "equal-weight mean", overlap_threshold = overlapThreshold,
        source_pixels = source$cell, weights = mapping$weights,
        target_metadata = S4Vectors::metadata(stData)
    )
    if (!annotations) {
        S4Vectors::metadata(output)$mz_annotation <- NULL
        S4Vectors::metadata(output)$db_3 <- NULL
    }
    if (dropUnmapped) output <- output[, output$mapped, drop = FALSE]
    verbose_message(paste(
        "Mapped", sum(output$mapped), "of", ncol(output),
        "target pixels."
    ), verbose = verbose)
    output
}
