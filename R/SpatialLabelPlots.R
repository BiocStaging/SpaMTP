# Development assistance: see inst/COPYRIGHTS.

#' Plot pixel labels at their spatial positions
#'
#' Draws one point per pixel at `SpatialExperiment::spatialCoords(x)`,
#' coloured by a categorical `colData()` column such as tissue regions or
#' clusters, and optionally shaped by a second column. The y axis is reversed
#' to follow image orientation; multiple samples are shown in separate panels.
#'
#' @param x A SpatialExperiment, or an aligned Cardinal `MSImagingExperiment`.
#' @param colourBy Name of the `colData()` column used for point colour.
#' @param shapeBy Optional name of a `colData()` column used for point shape.
#' @param colours Optional colours, named by label, or one per label in level
#'   order. NULL uses `grDevices::hcl.colors(n, "Dark 3")`.
#' @param pointSize Point size.
#' @param alpha Point opacity.
#' @param legendTitle Title of the colour legend (default `colourBy`).
#'
#' @return A `ggplot`.
#' @export
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(1, 1, 6, dimnames = list("mz-100", NULL))),
#'     colData = S4Vectors::DataFrame(
#'         region = rep(c("cortex", "medulla"), 3),
#'         side = rep(c("left", "right"), each = 3)
#'     ),
#'     spatialCoords = cbind(x = rep(1:3, 2), y = rep(1:2, each = 3))
#' )
#' plotSpatialLabels(spe, "region", shapeBy = "side", pointSize = 3)
plotSpatialLabels <- function(
    x, colourBy, shapeBy = NULL, colours = NULL,
    pointSize = 1, alpha = 1, legendTitle = colourBy
) {
    if (methods::is(x, "MSImagingExperiment")) {
        x <- asSpatialExperiment(x)
    }
    coordinates <- .nativeCoordinates(x, requireNames = FALSE)
    metadata <- SummarizedExperiment::colData(x)
    columns <- c(colourBy, shapeBy)
    if (length(colourBy) != 1L || length(shapeBy) > 1L ||
        !all(columns %in% names(metadata))) {
        stop(
            "colourBy and shapeBy must each name one colData() column.",
            call. = FALSE
        )
    }
    labels <- metadata[[colourBy]]
    if (!is.factor(labels)) {
        labels <- factor(
            labels,
            levels = .sortLabels(unique(stats::na.omit(labels)))
        )
    }
    coordinates$label <- labels
    if (is.null(colours)) {
        colours <- stats::setNames(
            grDevices::hcl.colors(nlevels(labels), "Dark 3"),
            levels(labels)
        )
    } else if (is.null(names(colours))) {
        if (length(colours) < nlevels(labels)) {
            stop("Supply at least one colour per label.", call. = FALSE)
        }
        colours <- stats::setNames(
            colours[seq_len(nlevels(labels))],
            levels(labels)
        )
    } else if (!all(levels(labels) %in% names(colours))) {
        stop("Named colours must cover every label.", call. = FALSE)
    }
    mapping <- ggplot2::aes(x = .data$x, y = .data$y, colour = .data$label)
    if (!is.null(shapeBy)) {
        coordinates$shape <- factor(metadata[[shapeBy]])
        mapping <- ggplot2::aes(
            x = .data$x, y = .data$y, colour = .data$label,
            shape = .data$shape
        )
    }
    plot <- ggplot2::ggplot(coordinates, mapping) +
        ggplot2::geom_point(size = pointSize, alpha = alpha) +
        ggplot2::coord_fixed() +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_colour_manual(
            values = colours,
            na.value = "grey80",
            drop = FALSE
        ) +
        ggplot2::guides(colour = ggplot2::guide_legend(
            override.aes = list(size = 3, alpha = 1)
        )) +
        ggplot2::labs(colour = legendTitle, shape = shapeBy) +
        .spatialPlotTheme() +
        ggplot2::theme(legend.box = "vertical")
    if (length(unique(coordinates$sample_id)) > 1L) {
        plot <- plot + ggplot2::facet_wrap(~sample_id)
    }
    plot
}
