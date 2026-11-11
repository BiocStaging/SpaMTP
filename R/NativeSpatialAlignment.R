.alignmentMatrix <- function(angle, scaleX, scaleY, shiftX, shiftY, center) {
    values <- c(angle, scaleX, scaleY, shiftX, shiftY, center)
    if (length(center) != 2L || any(!is.finite(values)) || scaleX == 0 || scaleY == 0) {
        stop(
            "Alignment parameters must be finite with non-zero scales.",
            call. = FALSE
        )
    }
    angle <- angle * pi / 180
    linear <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L) %*%
        diag(c(scaleX, scaleY))
    output <- diag(3L)
    output[1:2, 1:2] <- linear
    output[1:2, 3] <- center + c(shiftX, shiftY) - as.numeric(linear %*% center)
    output
}

.singleAlignmentSample <- function(object) {
    samples <- unique(.nativeCoordinates(object)$sample_id)
    if (length(samples) != 1L) {
        stop(
            "Align one tissue sample at a time by subsetting the object.",
            call. = FALSE
        )
    }
    samples[[1L]]
}

.alignmentApp <- function(
    source, target, imageId = NULL, multiplier = 1,
    showTarget = TRUE
) {
    sourceFrame <- .nativeCoordinates(source)
    targetFrame <- .nativeCoordinates(target)
    sampleId <- .singleAlignmentSample(target)
    .singleAlignmentSample(source)
    image <- if (is.null(imageId)) NULL else .nativeImage(target, imageId, sampleId)
    initial <- diag(c(multiplier, multiplier, 1))
    center <- colMeans(sourceFrame[, c("x", "y")]) * multiplier
    extent <- max(
        diff(range(c(sourceFrame$x * multiplier, targetFrame$x))),
        diff(range(c(sourceFrame$y * multiplier, targetFrame$y))), 1
    )
    ui <- shiny::fluidPage(
        shiny::titlePanel("Align spatial measurements"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::sliderInput("angle", "Rotation (degrees)", -180, 180, 0),
                shiny::sliderInput(
                    "shiftX",
                    "Horizontal shift",
                    -extent,
                    extent,
                    0
                ),
                shiny::sliderInput(
                    "shiftY",
                    "Vertical shift",
                    -extent,
                    extent,
                    0
                ),
                shiny::sliderInput(
                    "scaleX",
                    "Horizontal scale",
                    0.1,
                    4,
                    1,
                    step = 0.01
                ),
                shiny::sliderInput(
                    "scaleY",
                    "Vertical scale",
                    0.1,
                    4,
                    1,
                    step = 0.01
                ),
                shiny::checkboxInput("flipX", "Flip horizontally", FALSE),
                shiny::checkboxInput("flipY", "Flip vertically", FALSE),
                shiny::actionButton("apply", "Return aligned data"),
                shiny::actionButton("cancel", "Cancel")
            ),
            shiny::mainPanel(shiny::plotOutput("preview", height = "700px"))
        )
    )
    server <- function(input, output, session) {
        transform <- shiny::reactive({
            .alignmentMatrix(
                input$angle,
                input$scaleX * if (isTRUE(input$flipX)) -1 else 1,
                input$scaleY * if (isTRUE(input$flipY)) -1 else 1,
                input$shiftX, input$shiftY, center
            ) %*% initial
        })
        output$preview <- shiny::renderPlot({
            xy <- cbind(as.matrix(sourceFrame[, c("x", "y")]), 1) %*% t(transform())
            frame <- data.frame(x = xy[, 1], y = xy[, 2])
            plot <- ggplot2::ggplot()
            if (!is.null(image)) plot <- plot + .imageLayer(image)
            if (showTarget) {
                plot <- plot + ggplot2::geom_point(
                    data = targetFrame, ggplot2::aes(x = .data$x, y = .data$y), colour = "steelblue", alpha = 0.5
                )
            }
            plot + ggplot2::geom_point(
                data = frame, ggplot2::aes(x = .data$x, y = .data$y),
                colour = "tomato", alpha = 0.7
            ) + ggplot2::coord_fixed() +
                ggplot2::scale_y_reverse() + ggplot2::theme_void()
        })
        shiny::observeEvent(input$apply, shiny::stopApp(transform()))
        shiny::observeEvent(input$cancel, shiny::stopApp(NULL))
    }
    shiny::shinyApp(ui, server)
}

#' Align two spatial experiments using an affine transformation
#'
#' A supplied homogeneous matrix is applied reproducibly without a GUI. With
#' transformation = NULL, an interactive preview allows translation, rotation,
#' scaling and reflection. Coordinates are always stored in their full spatial
#' units; an imgData scaleFactor affects only display of the reference image.
#' @param smData,stData Source and reference SpatialExperiment objects, each
#'   containing one tissue sample.
#' @param transformation A 3 by 3 homogeneous matrix mapping source coordinates
#'   to reference coordinates, or NULL for interactive alignment.
#' @param imageSlice Optional reference image ID in `imgData(stData)`.
#' @param msiPixelMultiplier Explicit initial scaling of source coordinates.
#'   Defaults to one. Applied before transformation, including GUI transforms.
#' @param shinyHost,shinyPort Host and optional port for the interactive app.
#' @param verbose Show alignment messages.
#' @return The source SpatialExperiment with updated spatialCoords and
#'   alignment provenance recorded by applySpatialAlignment.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' transform <- diag(3)
#' transform[1, 3] <- 2
#' aligned <- alignSpatialOmics(spe, spe, transformation = transform)
#' checkAlignment(aligned, spe)
alignSpatialOmics <- function(
    smData, stData, transformation = NULL, imageSlice = NULL,
    msiPixelMultiplier = 1, shinyHost = "127.0.0.1", shinyPort = NULL,
    verbose = FALSE
) {
    smData <- .nativeSpatialObject(smData)
    stData <- .nativeSpatialObject(stData)
    .singleAlignmentSample(smData)
    .singleAlignmentSample(stData)
    if (length(msiPixelMultiplier) != 1L || !is.finite(msiPixelMultiplier) ||
        msiPixelMultiplier <= 0) {
        stop("msiPixelMultiplier must be one positive number.", call. = FALSE)
    }
    if (is.null(transformation)) {
        transformation <- shiny::runApp(
            .alignmentApp(smData, stData, imageSlice, msiPixelMultiplier),
            host = shinyHost, port = shinyPort
        )
        if (is.null(transformation)) stop("Alignment cancelled.", call. = FALSE)
    } else {
        if (!is.matrix(transformation) || !identical(dim(transformation), c(3L, 3L)) ||
            any(!is.finite(transformation))) {
            stop(
                "transformation must be a finite 3 by 3 matrix.",
                call. = FALSE
            )
        }
        transformation <- transformation %*% diag(c(
            msiPixelMultiplier,
            msiPixelMultiplier, 1
        ))
    }
    applySpatialAlignment(smData,
        stData = stData,
        alignment = list(transformation = list(matrix = transformation)),
        verbose = verbose
    )
}

#' Inspect the alignment of two spatial experiments
#' @param smData,stData SpatialExperiment objects in the same coordinate
#'   units.
#' @param names Two labels for the source and target measurements.
#' @param cols Two point colours.
#' @param imageSlice Optional image ID in ST.data. With multiple samples,
#'   subset to one sample before drawing an image overlay.
#' @param size Point size.
#' @return A ggplot; multiple samples are shown in separate panels.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' checkAlignment(spe, spe, names = c("source", "reference"))
checkAlignment <- function(
    smData, stData, names = c("SM", "ST"), cols = NULL,
    imageSlice = NULL, size = 0.5
) {
    source <- .nativeCoordinates(smData)
    target <- .nativeCoordinates(stData)
    if (length(names) != 2L || (!is.null(cols) && length(cols) != 2L)) {
        stop("Provide two dataset labels and two colours.", call. = FALSE)
    }
    source$dataset <- names[[1L]]
    target$dataset <- names[[2L]]
    frame <- rbind(source, target)
    plot <- ggplot2::ggplot(
        frame,
        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$dataset)
    )
    if (!is.null(imageSlice)) {
        .singleAlignmentSample(smData)
        sampleId <- .singleAlignmentSample(stData)
        plot <- plot + .imageLayer(.nativeImage(stData, imageSlice, sampleId))
    }
    plot <- plot + ggplot2::geom_point(size = size) +
        ggplot2::scale_colour_manual(values = cols %||% c("#F8766D", "#00BFC4")) +
        ggplot2::scale_y_reverse() + ggplot2::coord_fixed() + ggplot2::theme_void()
    if (is.null(imageSlice) && length(unique(frame$sample_id)) > 1L) {
        plot <- plot + ggplot2::facet_wrap(~sample_id)
    }
    plot
}

#' Attach and optionally align an optical image to MSI
#' @param imagePath Local image path.
#' @param object A SpatialExperiment.
#' @param imageId Identifier for the attached image.
#' @param sampleId Sample identifier; required for multi-sample objects.
#' @param scaleFactor Image pixels per spatial coordinate unit.
#' @param transformation Optional 3 by 3 matrix mapping MSI to image coordinates
#'   in full spatial units.
#' @param interactive Open the manual alignment app when no matrix is supplied.
#' @param ... Additional arguments to shiny::runApp when interactive is TRUE.
#' @return A SpatialExperiment with an imgData entry and optionally transformed
#'   spatialCoords. Coordinate alignment requires one tissue sample.
#' @export
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' image_file <- tempfile(fileext = ".png")
#' png::writePNG(matrix(seq(0, 1, length.out = 24), 4, 6), image_file)
#' spe <- addSMImage(image_file, spe, interactive = FALSE)
#' SpatialExperiment::imgData(spe)
#' unlink(image_file)
addSMImage <- function(
    imagePath, object, imageId = "optical", sampleId = NULL,
    scaleFactor = 1, transformation = NULL,
    interactive = TRUE, ...
) {
    object <- .nativeSpatialObject(object)
    sampleId <- sampleId %||% .singleAlignmentSample(object)
    reference <- addSpatialImage(object,
        imageSource = imagePath,
        sampleId = sampleId, imageId = imageId, scaleFactor = scaleFactor
    )
    if (interactive && is.null(transformation)) {
        transformation <- shiny::runApp(
            .alignmentApp(object, reference, imageId, showTarget = FALSE), ...
        )
        if (is.null(transformation)) stop("Alignment cancelled.", call. = FALSE)
    }
    if (!is.null(transformation)) {
        object <- alignSpatialOmics(
            object,
            reference,
            transformation = transformation
        )
        SpatialExperiment::imgData(object) <- SpatialExperiment::imgData(reference)
        return(object)
    }
    reference
}
