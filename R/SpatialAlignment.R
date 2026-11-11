#' Apply automatic or precomputed spatial alignment
#'
#' Aligns a moving Spatial Metabolomics (SM) `SpatialExperiment`
#' to a fixed SpatialExperiment. The function can apply coordinates exported by
#' the SMINT workflow, apply a homogeneous transformation matrix, estimate an
#' affine transform from paired landmarks, or run the STalign LDDMM backend used
#' by SMINT through `reticulate`.
#'
#' Python is only required when `method = "lddmm"` and `alignment = NULL`.
#' Precomputed SMINT coordinates and affine landmark alignment are handled
#' entirely in R. The backend summary records STalign's affine component as
#' `affine_yx`; it is not the complete nonlinear LDDMM transformation.
#'
#' @param smData A `SpatialExperiment` containing the
#'   moving SM coordinates.
#' @param stData An optional matching spatial container with fixed coordinates.
#'   It is required when computing an affine or LDDMM alignment, but is optional
#'   when applying final coordinates supplied through `alignment`.
#' @param alignment Optional precomputed alignment. Accepted inputs are a data
#'   frame or matrix of final coordinates; a 2-by-3 or 3-by-3 homogeneous
#'   transformation matrix; a list containing `aligned_coordinates`,
#'   `coordinates`, `transformation`, or `matrix`; a CSV/TSV/RDS/JSON file; or a
#'   SMINT output directory containing `aligned_coordinates.csv` and optionally
#'   `transformation.json`. Coordinate columns produced by the existing SpaMTP
#'   and SMINT notebooks (`x_transformed`, `x_new`, and `x_final`, with their y
#'   counterparts) are detected automatically.
#' @param method Alignment method used when `alignment` is `NULL`. `"affine"`
#'   estimates a six-degree-of-freedom transformation from paired landmarks in
#'   R. `"lddmm"` runs STalign's landmark-guided affine plus diffeomorphic
#'   registration through Python.
#' @param smFov Optional imgData image ID for the moving experiment. By
#'   default,
#'   the first image in `smData` is used.
#' @param stImage Optional imgData image ID for the fixed experiment. By
#'   default, the first image in `stData` is used.
#' @param smLandmarks,stLandmarks Matched landmark coordinates as two-column
#' matrices or data frames in x-y order. At least three non-collinear pairs are
#'   required for `method = "affine"`. Landmarks are optional but strongly
#'   recommended for `method = "lddmm"` when the sections are not already close.
#'   ST landmarks must use the target coordinate system after any
#'   `stScaleFactor` has been applied.
#' @param landmarkOrder Coordinate-column order in both landmark inputs.
#'   `"xy"` is the R/SpaMTP convention. Use `"yx"` for row-column arrays saved
#'   by the notebook's point annotator.
#' @param smLandmarkSpace Whether SM landmarks use the `"original"` FOV
#'   coordinates or the `"preprocessed"` coordinates after source scaling,
#'   rotation, reflection, translation, and origin adjustment. Point-annotator
#'   landmarks made from the rasterized notebook output are preprocessed.
#' @param coordinateColumns Optional names or integer positions of the x and y
#'   columns in a precomputed coordinate table.
#' @param cellColumn Optional name or integer position of the cell/spot ID
#'   column in a precomputed coordinate table. When no IDs are available, rows
#'   must already be in the same order as the moving SM centroids.
#' @param stScaleFactor Optional fixed-coordinate scale factor. Supply a
#'   numeric value or a scale-factor name such as `"hires"` or `"lowres"` for
#'   Visium images. `NULL` uses full-resolution coordinates.
#' @param sourceScale Numeric scalar or x-y vector applied to moving
#'   coordinates
#'   before fitting. This reproduces the initial SM scaling step in the SMINT
#'   notebooks.
#' @param sourceRotate Counter-clockwise rotation in degrees, around the
#'   origin,
#'   applied after scaling.
#' @param sourceFlip Logical x-y vector indicating reflection around the
#'   origin.
#' @param sourceTranslate Numeric x-y translation applied after scaling,
#'   reflection, and rotation.
#' @param sourceOrigin Optional numeric x-y lower bound. Coordinates are
#'   shifted
#'   only where necessary so that their minima are at least this value. Use
#'   `c(0, 0)` to reproduce the non-negative coordinate step in the notebook.
#' @param dx Raster spacing passed to `STalign.rasterize()`.
#' @param niter Number of LDDMM optimisation iterations.
#' @param epV Velocity-field gradient step passed to `STalign.LDDMM()`.
#' @param device One of `"auto"`, `"cpu"`, or `"cuda"`. `"auto"` uses CUDA
#'   when Torch reports that it is available.
#' @param seed Random seed set for NumPy and Torch before LDDMM.
#' @param lddmmArgs Named list of additional arguments passed to
#'   `STalign.LDDMM()`. Values here override `niter` and `epV`.
#' @param python Optional path to a Python executable containing `STalign` and
#'   `torch`. It must be selected before reticulate initialises Python.
#' @param alignmentName Name used to store alignment provenance in
#'   `S4Vectors::metadata(smData)$spatial_alignment`.
#' @param store Logical; store lightweight provenance and diagnostics in the
#'   returned object.
#' @param return Either `"object"` (default) or `"result"`. The latter returns a
#'   list containing the object, coordinate tables, diagnostics, parameters, and
#'   lightweight backend information.
#' @param diagnosticsMaxPoints Maximum number of moving and fixed points used
#'   for nearest-neighbour diagnostics.
#' @param verbose Logical; print progress messages.
#'
#' @return The updated spatial container when `return = "object"`; otherwise an
#'   object of class
#'   `spamtp_spatial_alignment` containing the aligned object and alignment
#'   details.
#'
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' transform <- diag(3)
#' transform[1:2, 3] <- c(10, 20)
#' aligned <- applySpatialAlignment(spe,
#'     alignment = list(transformation = list(matrix = transform)), verbose = FALSE
#' )
#' head(SpatialExperiment::spatialCoords(aligned))
applySpatialAlignment <- function(
    smData,
    stData = NULL,
    alignment = NULL,
    method = c("lddmm", "affine"),
    smFov = NULL,
    stImage = NULL,
    smLandmarks = NULL,
    stLandmarks = NULL,
    landmarkOrder = c("xy", "yx"),
    smLandmarkSpace = c("original", "preprocessed"),
    coordinateColumns = NULL,
    cellColumn = NULL,
    stScaleFactor = NULL,
    sourceScale = 1,
    sourceRotate = 0,
    sourceFlip = c(FALSE, FALSE),
    sourceTranslate = c(0, 0),
    sourceOrigin = NULL,
    dx = 30,
    niter = 1000,
    epV = 200,
    device = c("auto", "cpu", "cuda"),
    seed = 1,
    lddmmArgs = list(),
    python = NULL,
    alignmentName = "SMINT",
    store = TRUE,
    return = c("object", "result"),
    diagnosticsMaxPoints = 10000,
    verbose = TRUE
) {
    method <- match.arg(method)
    device <- match.arg(device)
    return <- match.arg(return)
    landmarkOrder <- match.arg(landmarkOrder)
    smLandmarkSpace <- match.arg(smLandmarkSpace)
    if (!is.character(alignmentName) || length(alignmentName) != 1L ||
        is.na(alignmentName) || !nzchar(alignmentName)) {
        stop(
            "`alignmentName` must be one non-empty character string.",
            call. = FALSE
        )
    }
    for (argument in c("store", "verbose")) {
        value <- get(argument)
        if (!is.logical(value) || length(value) != 1L || is.na(value)) {
            stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
        }
    }

    .sa_validate_spatial_container(smData, "smData")
    smFov <- .sa_resolve_image(smData, smFov, "smFov")
    source <- .sa_get_coordinates(smData, smFov)

    preprocessing <- .sa_preprocess_coordinates(
        source[, c("x", "y"), drop = FALSE],
        scale = sourceScale,
        rotate = sourceRotate,
        flip = sourceFlip,
        translate = sourceTranslate,
        origin = sourceOrigin
    )
    moving <- preprocessing$coordinates
    moving$cell <- source$cell

    target <- NULL
    if (!is.null(stData)) {
        .sa_validate_spatial_container(stData, "stData")
        stImage <- .sa_resolve_image(stData, stImage, "stImage")
        target <- .sa_get_coordinates(stData, stImage)
        target_scale <- .sa_scale_factor(stData, stImage, stScaleFactor)
        target[, c("x", "y")] <- target[, c("x", "y"), drop = FALSE] * target_scale
    }

    landmarks <- .sa_prepare_landmarks(
        smLandmarks = smLandmarks,
        stLandmarks = stLandmarks,
        preprocessing = preprocessing$matrix,
        landmarkOrder = landmarkOrder,
        smLandmarkSpace = smLandmarkSpace
    )

    backend <- "precomputed"
    backend_info <- list()
    transformation <- NULL
    transformed_landmarks <- NULL

    if (!is.null(alignment)) {
        payload <- .sa_read_alignment(alignment)
        coordinates_payload <- .sa_find_coordinates(payload)
        transformation <- .sa_find_transformation(payload)

        if (!is.null(coordinates_payload)) {
            aligned <- .sa_match_coordinates(
                coordinates_payload,
                sourceCells = source$cell,
                coordinateColumns = coordinateColumns,
                cellColumn = cellColumn
            )
        } else if (!is.null(transformation)) {
            aligned_xy <- .sa_apply_transform(
                moving[, c("x", "y"), drop = FALSE],
                transformation
            )
            aligned <- data.frame(
                cell = source$cell,
                x = aligned_xy[, 1],
                y = aligned_xy[, 2]
            )
            if (!is.null(landmarks)) {
                transformed_landmarks <- .sa_apply_transform(landmarks$source, transformation)
            }
        } else {
            stop(
                "`alignment` did not contain aligned coordinates or a supported transformation matrix.",
                call. = FALSE
            )
        }
    } else {
        if (is.null(target)) {
            stop(
                "`stData` is required when `alignment` is NULL.",
                call. = FALSE
            )
        }

        if (identical(method, "affine")) {
            if (is.null(landmarks)) {
                stop(
                    "Paired `smLandmarks` and `stLandmarks` are required for affine alignment.",
                    call. = FALSE
                )
            }
            transformation <- .sa_fit_affine(landmarks$source, landmarks$target)
            aligned_xy <- .sa_apply_transform(
                moving[, c("x", "y"), drop = FALSE],
                transformation
            )
            transformed_landmarks <- .sa_apply_transform(landmarks$source, transformation)
            aligned <- data.frame(
                cell = source$cell,
                x = aligned_xy[, 1],
                y = aligned_xy[, 2]
            )
            backend <- "R-affine"
        } else {
            if (verbose) {
                message("Running SMINT-compatible STalign LDDMM registration ...")
            }
            fit <- .sa_run_stalign(
                moving = moving[, c("x", "y"), drop = FALSE],
                target = target[, c("x", "y"), drop = FALSE],
                landmarks = landmarks,
                dx = dx,
                niter = niter,
                epV = epV,
                device = device,
                seed = seed,
                lddmmArgs = lddmmArgs,
                python = python
            )
            aligned <- data.frame(
                cell = source$cell,
                x = fit$coordinates[, 1],
                y = fit$coordinates[, 2]
            )
            transformation <- NULL
            transformed_landmarks <- fit$landmarks
            backend_info <- fit$backend
            backend <- "STalign-LDDMM"
        }
    }

    if (any(!is.finite(as.matrix(aligned[, c("x", "y"), drop = FALSE])))) {
        stop("Aligned coordinates contain non-finite values.", call. = FALSE)
    }

    diagnostics <- .sa_alignment_diagnostics(
        original = source[, c("x", "y"), drop = FALSE],
        aligned = aligned[, c("x", "y"), drop = FALSE],
        target = if (is.null(target)) NULL else target[, c("x", "y"), drop = FALSE],
        transformedLandmarks = transformed_landmarks,
        targetLandmarks = if (is.null(landmarks)) NULL else landmarks$target,
        maxPoints = diagnosticsMaxPoints
    )

    output <- .sa_set_fov_coordinates(
        smData,
        aligned = aligned
    )

    parameters <- list(
        method = if (is.null(alignment)) method else "precomputed",
        backend = backend,
        SM.fov = smFov,
        ST.image = stImage,
        ST.scale.factor = stScaleFactor,
        source.scale = sourceScale,
        source.rotate = sourceRotate,
        source.flip = sourceFlip,
        source.translate = sourceTranslate,
        source.origin = sourceOrigin,
        landmark.order = landmarkOrder,
        SM.landmark.space = smLandmarkSpace,
        dx = if (identical(backend, "STalign-LDDMM")) dx else NULL,
        niter = if (identical(backend, "STalign-LDDMM")) niter else NULL,
        epV = if (identical(backend, "STalign-LDDMM")) epV else NULL,
        device = if (identical(backend, "STalign-LDDMM")) device else NULL,
        seed = if (identical(backend, "STalign-LDDMM")) seed else NULL,
        lddmm.args = if (identical(backend, "STalign-LDDMM")) lddmmArgs else NULL
    )

    provenance <- list(
        name = alignmentName,
        engine = "SpaMTP::applySpatialAlignment",
        backend = backend,
        parameters = parameters,
        preprocessing_matrix = preprocessing$matrix,
        transformation = transformation,
        diagnostics = diagnostics,
        backend_info = backend_info
    )

    if (isTRUE(store)) {
        storedAlignments <- .storedData(output, "spatial_alignment") %||% list()
        storedAlignments[[alignmentName]] <- provenance
        output <- .setStoredData(output, "spatial_alignment", storedAlignments)
    }

    coordinates <- data.frame(
        cell = source$cell,
        x_original = source$x,
        y_original = source$y,
        x_aligned = aligned$x,
        y_aligned = aligned$y,
        stringsAsFactors = FALSE
    )

    result <- structure(
        list(
            object = output,
            coordinates = coordinates,
            target_coordinates = target,
            diagnostics = diagnostics,
            parameters = parameters,
            preprocessing_matrix = preprocessing$matrix,
            transformation = transformation,
            backend = backend_info,
            provenance = provenance
        ),
        class = "spamtp_spatial_alignment"
    )

    if (verbose) {
        detail <- if (!is.null(diagnostics$nearest_target_after$median)) {
            paste0(
                "; median nearest-target distance ",
                signif(diagnostics$nearest_target_before$median, 4), " -> ",
                signif(diagnostics$nearest_target_after$median, 4)
            )
        } else {
            ""
        }
        message(
            "Aligned ", nrow(aligned), " SM coordinates using ", backend,
            detail
        )
    }

    if (identical(return, "object")) output else result
}


.sa_validate_spatial_container <- function(object, argument) {
    .requireExperiment(object, "SpatialExperiment")
}


.sa_resolve_image <- function(object, image, argument) {
    .requireExperiment(object, "SpatialExperiment")
    images <- as.character(SpatialExperiment::imgData(object)$image_id)
    if (is.null(image)) {
        return(if (length(images)) images[[1L]] else NULL)
    }
    if (!length(images) || !image %in% images) {
        stop(
            "Invalid `", argument, "`. Available imgData values: ",
            paste(images, collapse = ", "), ".",
            call. = FALSE
        )
    }
    return(image)
}


.sa_get_coordinates <- function(object, image) {
    coordinates <- as.data.frame(SpatialExperiment::spatialCoords(object))
    coordinates$cell <- colnames(object)
    if (!all(c("x", "y") %in% names(coordinates))) {
        stop(
            "Spatial coordinates for `",
            image,
            "` do not contain x and y columns.",
            call. = FALSE
        )
    }
    if (!"cell" %in% names(coordinates)) {
        cells <- rownames(coordinates)
        if (is.null(cells) || any(!nzchar(cells))) {
            stop(
                "Spatial coordinates for `",
                image,
                "` do not contain cell identifiers.",
                call. = FALSE
            )
        }
        coordinates$cell <- cells
    }
    coordinates <- coordinates[, c("x", "y", "cell"), drop = FALSE]
    coordinates$x <- suppressWarnings(as.numeric(coordinates$x))
    coordinates$y <- suppressWarnings(as.numeric(coordinates$y))
    coordinates$cell <- as.character(coordinates$cell)
    if (any(!is.finite(as.matrix(coordinates[, c("x", "y"), drop = FALSE])))) {
        stop(
            "Spatial coordinates for `",
            image,
            "` contain non-finite values.",
            call. = FALSE
        )
    }
    if (anyDuplicated(coordinates$cell)) {
        stop(
            "Spatial coordinates for `",
            image,
            "` contain duplicated cell identifiers.",
            call. = FALSE
        )
    }
    rownames(coordinates) <- NULL
    coordinates
}


.sa_scale_factor <- function(object, image, scaleFactor) {
    if (is.null(scaleFactor)) {
        return(1)
    }
    if (is.numeric(scaleFactor)) {
        if (length(scaleFactor) != 1L || !is.finite(scaleFactor) || scaleFactor <= 0) {
            stop(
                "Numeric `stScaleFactor` must be one finite value greater than zero.",
                call. = FALSE
            )
        }
        return(as.numeric(scaleFactor))
    }
    if (!is.character(scaleFactor) || length(scaleFactor) != 1L || is.na(scaleFactor)) {
        stop(
            "`stScaleFactor` must be NULL, numeric, or one scale-factor name.",
            call. = FALSE
        )
    }
    imageData <- as.data.frame(SpatialExperiment::imgData(object))
    selected <- which(imageData$image_id == image)
    if (length(selected) != 1L) {
        stop(
            "A unique imgData row is required for the selected image.",
            call. = FALSE
        )
    }
    if (!scaleFactor %in% c("scaleFactor", "hires", "lowres")) {
        stop(
            "SpatialExperiment images provide a single `scaleFactor`; use that name, `hires`, or `lowres`.",
            call. = FALSE
        )
    }
    value <- as.numeric(imageData$scaleFactor[[selected]])
    if (!is.finite(value) || value <= 0) {
        stop(
            "Selected image scale factor must be finite and greater than zero.",
            call. = FALSE
        )
    }
    return(value)
}


.sa_pair <- function(value, argument, mode = c("numeric", "logical")) {
    mode <- match.arg(mode)
    if (length(value) == 1L) value <- rep(value, 2L)
    if (length(value) != 2L || anyNA(value)) {
        stop(
            "`",
            argument,
            "` must contain one or two non-missing values.",
            call. = FALSE
        )
    }
    if (identical(mode, "numeric")) {
        value <- suppressWarnings(as.numeric(value))
        if (any(!is.finite(value))) {
            stop(
                "`",
                argument,
                "` must contain finite numeric values.",
                call. = FALSE
            )
        }
    } else {
        if (!is.logical(value)) {
            stop("`", argument, "` must be logical.", call. = FALSE)
        }
    }
    value
}


.sa_apply_transform <- function(points, transformation) {
    points <- as.matrix(points)
    storage.mode(points) <- "double"
    if (ncol(points) != 2L || any(!is.finite(points))) {
        stop(
            "Points must be a finite numeric matrix with two columns.",
            call. = FALSE
        )
    }
    transformation <- .sa_matrix(transformation)
    if (identical(dim(transformation), c(2L, 3L))) {
        transformation <- rbind(transformation, c(0, 0, 1))
    }
    if (!identical(dim(transformation), c(3L, 3L))) {
        stop(
            "Transformation must be a 2-by-3 or 3-by-3 numeric matrix.",
            call. = FALSE
        )
    }
    if (any(!is.finite(transformation))) {
        stop(
            "Transformation must contain only finite numeric values.",
            call. = FALSE
        )
    }
    homogeneous <- cbind(points, 1)
    transformed <- homogeneous %*% t(transformation)
    denominator <- transformed[, 3L]
    if (any(!is.finite(denominator)) || any(abs(denominator) < .Machine$double.eps)) {
        stop(
            "Transformation produced invalid homogeneous coordinates.",
            call. = FALSE
        )
    }
    transformed[, c(1L, 2L), drop = FALSE] / denominator
}


.sa_preprocess_coordinates <- function(
    coordinates,
    scale = 1,
    rotate = 0,
    flip = c(FALSE, FALSE),
    translate = c(0, 0),
    origin = NULL
) {
    scale <- .sa_pair(scale, "sourceScale")
    if (any(scale == 0)) stop("`sourceScale` cannot contain zero.", call. = FALSE)
    flip <- .sa_pair(flip, "sourceFlip", mode = "logical")
    translate <- .sa_pair(translate, "sourceTranslate")
    if (!is.numeric(rotate) || length(rotate) != 1L || !is.finite(rotate)) {
        stop("`sourceRotate` must be one finite number.", call. = FALSE)
    }

    angle <- rotate * pi / 180
    scaling <- diag(c(scale * ifelse(flip, -1, 1), 1))
    rotation <- matrix(
        c(cos(angle), sin(angle), 0, -sin(angle), cos(angle), 0, 0, 0, 1),
        nrow = 3L
    )
    translation <- diag(3L)
    translation[c(1L, 2L), 3L] <- translate
    transformation <- translation %*% rotation %*% scaling
    transformed <- .sa_apply_transform(coordinates, transformation)

    if (!is.null(origin)) {
        origin <- .sa_pair(origin, "sourceOrigin")
        shift <- pmax(origin - apply(transformed, 2L, min), 0)
        origin_transform <- diag(3L)
        origin_transform[c(1L, 2L), 3L] <- shift
        transformation <- origin_transform %*% transformation
        transformed <- .sa_apply_transform(coordinates, transformation)
    }

    colnames(transformed) <- c("x", "y")
    list(coordinates = as.data.frame(transformed), matrix = transformation)
}


.sa_prepare_landmarks <- function(
    smLandmarks,
    stLandmarks,
    preprocessing,
    landmarkOrder,
    smLandmarkSpace
) {
    if (is.null(smLandmarks) && is.null(stLandmarks)) {
        return(NULL)
    }
    if (is.null(smLandmarks) || is.null(stLandmarks)) {
        stop(
            "`smLandmarks` and `stLandmarks` must be supplied together.",
            call. = FALSE
        )
    }
    source <- .sa_landmark_matrix(smLandmarks, "smLandmarks")
    target <- .sa_landmark_matrix(stLandmarks, "stLandmarks")
    if (nrow(source) != nrow(target)) {
        stop(
            "`smLandmarks` and `stLandmarks` must contain the same number of rows.",
            call. = FALSE
        )
    }
    if (identical(landmarkOrder, "yx")) {
        source <- source[, c(2L, 1L), drop = FALSE]
        target <- target[, c(2L, 1L), drop = FALSE]
    }
    if (identical(smLandmarkSpace, "original")) {
        source <- .sa_apply_transform(source, preprocessing)
    }
    colnames(source) <- colnames(target) <- c("x", "y")
    list(source = source, target = target)
}


.sa_landmark_matrix <- function(x, argument) {
    if (inherits(x, "python.builtin.object")) {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop(
                "Package `reticulate` is required to convert Python landmarks.",
                call. = FALSE
            )
        }
        x <- reticulate::py_to_r(x)
    }
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    if (ncol(x) != 2L || nrow(x) < 1L || any(!is.finite(x))) {
        stop(
            "`",
            argument,
            "` must be a finite numeric matrix with two columns.",
            call. = FALSE
        )
    }
    colnames(x) <- c("x", "y")
    x
}


.sa_fit_affine <- function(source, target) {
    if (nrow(source) < 3L) {
        stop(
            "Affine alignment requires at least three paired landmarks.",
            call. = FALSE
        )
    }
    design <- cbind(source, 1)
    if (qr(design)$rank < 3L) {
        stop(
            "Affine landmarks must contain at least three non-collinear source points.",
            call. = FALSE
        )
    }
    coefficients <- qr.solve(design, target)
    rbind(t(coefficients), c(0, 0, 1))
}


.sa_read_alignment <- function(alignment) {
    if (inherits(alignment, "python.builtin.object")) {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop(
                "Package `reticulate` is required to convert a Python alignment object.",
                call. = FALSE
            )
        }
        alignment <- reticulate::py_to_r(alignment)
    }
    if (!is.character(alignment) || length(alignment) != 1L) {
        return(alignment)
    }
    if (!file.exists(alignment)) {
        stop("Alignment path does not exist: ", alignment, call. = FALSE)
    }
    if (dir.exists(alignment)) {
        coordinate_file <- file.path(alignment, "aligned_coordinates.csv")
        transform_file <- file.path(alignment, "transformation.json")
        if (!file.exists(coordinate_file) && !file.exists(transform_file)) {
            stop(
                "SMINT alignment directory must contain aligned_coordinates.csv or transformation.json.",
                call. = FALSE
            )
        }
        return(list(
            aligned_coordinates = if (file.exists(coordinate_file)) {
                utils::read.csv(coordinate_file, check.names = FALSE)
            } else {
                NULL
            },
            transformation = if (file.exists(transform_file)) {
                jsonlite::fromJSON(transform_file, simplifyVector = TRUE)
            } else {
                NULL
            }
        ))
    }

    extension <- tolower(sub("^.*[.]", "", basename(alignment)))
    switch(extension,
        csv = utils::read.csv(alignment, check.names = FALSE),
        tsv = utils::read.delim(alignment, check.names = FALSE),
        txt = utils::read.delim(alignment, check.names = FALSE),
        rds = readRDS(alignment),
        json = jsonlite::fromJSON(alignment, simplifyVector = TRUE),
        stop(
            "Unsupported alignment file extension: .",
            extension,
            call. = FALSE
        )
    )
}


.sa_find_coordinates <- function(payload) {
    if (is.data.frame(payload) || (is.matrix(payload) && !.sa_is_transform(payload))) {
        return(payload)
    }
    if (!is.list(payload)) {
        return(NULL)
    }
    for (name in c(
        "aligned_coordinates",
        "coordinates",
        "aligned.coords",
        "aligned"
    )) {
        if (!is.null(payload[[name]])) {
            return(payload[[name]])
        }
    }
    NULL
}


.sa_find_transformation <- function(payload) {
    if (.sa_is_transform(payload)) {
        return(.sa_matrix(payload))
    }
    if (!is.list(payload)) {
        return(NULL)
    }
    candidates <- list(
        payload$transformation,
        payload$transform,
        payload$matrix
    )
    for (candidate in candidates) {
        if (is.null(candidate)) next
        if (is.list(candidate) && !is.null(candidate$matrix)) candidate <- candidate$matrix
        if (.sa_is_transform(candidate)) {
            return(.sa_matrix(candidate))
        }
    }
    NULL
}


.sa_matrix <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
        if (!is.null(x$matrix)) x <- x$matrix
        if (is.list(x) && length(x) && all(vapply(x, is.atomic, logical(1)))) {
            x <- do.call(rbind, x)
        }
    }
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
}


.sa_is_transform <- function(x) {
    if (is.data.frame(x)) {
        return(FALSE)
    }
    candidate <- suppressWarnings(try(.sa_matrix(x), silent = TRUE))
    !inherits(candidate, "try-error") &&
        (identical(dim(candidate), c(2L, 3L)) || identical(dim(candidate), c(3L, 3L))) &&
        all(is.finite(candidate))
}


.sa_match_coordinates <- function(
    coordinates,
    sourceCells,
    coordinateColumns = NULL,
    cellColumn = NULL
) {
    coordinates <- as.data.frame(coordinates, check.names = FALSE)
    if (!nrow(coordinates)) stop(
        "Precomputed alignment contains no coordinates.",
        call. = FALSE
    )

    if (is.null(coordinateColumns)) {
        pairs <- list(
            c("x_final", "y_final"),
            c("x_transformed", "y_transformed"),
            c("x_new", "y_new"),
            c("x_aligned", "y_aligned"),
            c("aligned_x", "aligned_y"),
            c("x", "y")
        )
        coordinateColumns <- NULL
        for (pair in pairs) {
            if (all(pair %in% names(coordinates))) {
                coordinateColumns <- pair
                break
            }
        }
        if (is.null(coordinateColumns) && ncol(coordinates) == 2L) {
            coordinateColumns <- c(1L, 2L)
        }
    }
    if (length(coordinateColumns) != 2L ||
        (is.character(coordinateColumns) && !all(coordinateColumns %in% names(coordinates)))) {
        stop(
            "Could not identify aligned x-y columns. Supply `coordinateColumns = c(x, y)`.",
            call. = FALSE
        )
    }
    if (is.numeric(coordinateColumns) &&
        any(coordinateColumns < 1 | coordinateColumns > ncol(coordinates))) {
        stop(
            "Numeric `coordinateColumns` are outside the coordinate table.",
            call. = FALSE
        )
    }

    xy <- coordinates[, coordinateColumns, drop = FALSE]
    xy[] <- lapply(xy, function(value) suppressWarnings(as.numeric(value)))
    if (any(!is.finite(as.matrix(xy)))) {
        stop(
            "Precomputed aligned coordinates contain non-finite values.",
            call. = FALSE
        )
    }

    if (is.null(cellColumn)) {
        candidates <- c("cell", "barcode", "spot", "id", "Cell", "Barcode")
        matches <- candidates[candidates %in% names(coordinates)]
        cellColumn <- if (length(matches)) matches[[1L]] else NULL
    }

    ids <- NULL
    if (!is.null(cellColumn)) {
        if (length(cellColumn) != 1L ||
            (is.character(cellColumn) && !cellColumn %in% names(coordinates))) {
            stop("Invalid `cellColumn`.", call. = FALSE)
        }
        if (is.numeric(cellColumn) &&
            (cellColumn < 1 || cellColumn > ncol(coordinates))) {
            stop(
                "Numeric `cellColumn` is outside the coordinate table.",
                call. = FALSE
            )
        }
        ids <- as.character(coordinates[[cellColumn]])
    } else {
        candidate_rownames <- rownames(coordinates)
        default_rownames <- identical(
            candidate_rownames,
            as.character(seq_len(nrow(coordinates)))
        )
        if (!is.null(candidate_rownames) && !default_rownames && all(sourceCells %in% candidate_rownames)) {
            ids <- candidate_rownames
        }
    }

    if (is.null(ids)) {
        if (nrow(xy) != length(sourceCells)) {
            stop(
                "Precomputed coordinates have no cell IDs and their row count does not match the SM FOV.",
                call. = FALSE
            )
        }
        ids <- sourceCells
    }
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop(
            "Precomputed alignment contains missing or duplicated cell identifiers.",
            call. = FALSE
        )
    }
    index <- match(sourceCells, ids)
    if (anyNA(index)) {
        stop(
            "Precomputed alignment is missing ", sum(is.na(index)),
            " SM cell identifier(s).",
            call. = FALSE
        )
    }
    data.frame(
        cell = sourceCells,
        x = xy[[1L]][index],
        y = xy[[2L]][index],
        stringsAsFactors = FALSE
    )
}


.sa_set_fov_coordinates <- function(object, aligned) {
    .requireExperiment(object, "SpatialExperiment")
    index <- match(colnames(object), aligned$cell)
    if (anyNA(index)) {
        stop(
            "Aligned coordinates are missing pixels in the SpatialExperiment.",
            call. = FALSE
        )
    }
    coordinates <- SpatialExperiment::spatialCoords(object)
    coordinates[, c("x", "y")] <-
        as.matrix(aligned[index, c("x", "y"), drop = FALSE])
    SpatialExperiment::spatialCoords(object) <- coordinates
    return(object)
}


.sa_run_stalign <- function(
    moving,
    target,
    landmarks,
    dx,
    niter,
    epV,
    device,
    seed,
    lddmmArgs,
    python
) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package `reticulate` is required for `method = \"lddmm\"`. ",
            "Install it or apply precomputed SMINT coordinates.",
            call. = FALSE
        )
    }
    for (entry in list(dx = dx, niter = niter, epV = epV, seed = seed)) {
        if (!is.numeric(entry) || length(entry) != 1L || !is.finite(entry)) {
            stop(
                "`dx`, `niter`, `epV`, and `seed` must be finite scalar values.",
                call. = FALSE
            )
        }
    }
    if (dx <= 0 || niter < 1 || epV <= 0) {
        stop(
            "`dx`, `niter`, and `epV` must be greater than zero.",
            call. = FALSE
        )
    }
    if (niter != floor(niter) || seed != floor(seed)) {
        stop("`niter` and `seed` must be whole numbers.", call. = FALSE)
    }
    if (!is.list(lddmmArgs) || is.null(names(lddmmArgs)) && length(lddmmArgs)) {
        stop("`lddmmArgs` must be a named list.", call. = FALSE)
    }
    if (length(lddmmArgs) && any(!nzchar(names(lddmmArgs)))) {
        stop("Every entry in `lddmmArgs` must have a name.", call. = FALSE)
    }
    protected <- c(
        "xI",
        "I",
        "xJ",
        "J",
        "pointsI",
        "pointsJ",
        "L",
        "T",
        "A",
        "v",
        "xv",
        "device"
    )
    if (any(names(lddmmArgs) %in% protected)) {
        stop(
            "`lddmmArgs` cannot override internal arguments: ",
            paste(
                intersect(names(lddmmArgs), protected),
                collapse = ", "
            ), ".",
            call. = FALSE
        )
    }

    if (!is.null(python)) {
        if (!is.character(python) || length(python) != 1L || is.na(python) ||
            !file.exists(python)) {
            stop(
                "`python` must point to an existing Python executable.",
                call. = FALSE
            )
        }
        reticulate::use_python(python, required = TRUE)
    }
    if (!reticulate::py_module_available("STalign")) {
        stop(
            "Python module `STalign` is unavailable. Install STalign in the Python ",
            "environment selected by reticulate, or supply precomputed coordinates.",
            call. = FALSE
        )
    }
    if (!reticulate::py_module_available("torch")) {
        stop("Python module `torch` is required by STalign.", call. = FALSE)
    }

    stalign_module <- reticulate::import("STalign", convert = FALSE)
    stalign <- stalign_module$STalign
    np <- reticulate::import("numpy", convert = FALSE)
    torch <- reticulate::import("torch", convert = FALSE)

    np$random$seed(as.integer(seed))
    torch$manual_seed(as.integer(seed))
    cuda_available <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
    selected_device <- switch(device,
        auto = if (cuda_available) "cuda:0" else "cpu",
        cpu = "cpu",
        cuda = "cuda:0"
    )
    if (identical(device, "cuda") && !cuda_available) {
        stop(
            "`device = \"cuda\"` was requested but Torch reports no available CUDA device.",
            call. = FALSE
        )
    }
    if (cuda_available) torch$cuda$manual_seed_all(as.integer(seed))

    moving_x <- np$array(moving$x, dtype = "float64")
    moving_y <- np$array(moving$y, dtype = "float64")
    target_x <- np$array(target$x, dtype = "float64")
    target_y <- np$array(target$y, dtype = "float64")

    moving_raster <- stalign$rasterize(moving_x, moving_y, dx = as.numeric(dx))
    target_raster <- stalign$rasterize(target_x, target_y, dx = as.numeric(dx))
    xI <- reticulate::py_get_item(moving_raster, 0L)
    yI <- reticulate::py_get_item(moving_raster, 1L)
    I <- reticulate::py_get_item(moving_raster, 2L)
    xJ <- reticulate::py_get_item(target_raster, 0L)
    yJ <- reticulate::py_get_item(target_raster, 1L)
    J <- reticulate::py_get_item(target_raster, 2L)

    moving_yx <- np$stack(reticulate::tuple(moving_y, moving_x), axis = 1L)
    pointsI <- NULL
    pointsJ <- NULL
    initial_L <- NULL
    initial_T <- NULL
    if (!is.null(landmarks)) {
        pointsI <- np$array(
            landmarks$source[, c(2L, 1L), drop = FALSE],
            dtype = "float64"
        )
        pointsJ <- np$array(
            landmarks$target[, c(2L, 1L), drop = FALSE],
            dtype = "float64"
        )
        initial <- stalign$L_T_from_points(pointsI, pointsJ)
        initial_L <- reticulate::py_get_item(initial, 0L)
        initial_T <- reticulate::py_get_item(initial, 1L)
    }

    arguments <- list(
        reticulate::tuple(yI, xI),
        I,
        reticulate::tuple(yJ, xJ),
        J,
        pointsI = pointsI,
        pointsJ = pointsJ,
        L = initial_L,
        "T" = initial_T,
        niter = as.integer(niter),
        epV = as.numeric(epV),
        device = selected_device
    )
    arguments[names(lddmmArgs)] <- lddmmArgs

    fit <- tryCatch(
        do.call(stalign$LDDMM, arguments),
        error = function(error) {
            stop(
                "STalign LDDMM failed: ", conditionMessage(error),
                ". Check coordinate scale, landmarks, `dx`, and `lddmmArgs`.",
                call. = FALSE
            )
        }
    )

    transformed <- stalign$transform_points_source_to_target(
        fit$xv, fit$v, fit$A, moving_yx
    )
    transformed <- reticulate::py_to_r(transformed$detach()$cpu()$numpy())
    aligned <- transformed[, c(2L, 1L), drop = FALSE]
    colnames(aligned) <- c("x", "y")

    transformed_landmarks <- NULL
    if (!is.null(pointsI)) {
        transformed_landmarks <- stalign$transform_points_source_to_target(
            fit$xv, fit$v, fit$A, pointsI
        )
        transformed_landmarks <- reticulate::py_to_r(
            transformed_landmarks$detach()$cpu()$numpy()
        )[, c(2L, 1L), drop = FALSE]
        colnames(transformed_landmarks) <- c("x", "y")
    }

    affine <- reticulate::py_to_r(fit$A$detach()$cpu()$numpy())
    velocity_shape <- unlist(reticulate::py_to_r(fit$v$shape), use.names = FALSE)
    python_version <- reticulate::py_config()$version

    try(
        {
            pyplot <- reticulate::import("matplotlib.pyplot", convert = FALSE)
            pyplot$close("all")
        },
        silent = TRUE
    )

    list(
        coordinates = aligned,
        landmarks = transformed_landmarks,
        backend = list(
            engine = "STalign",
            python = as.character(python_version),
            device = selected_device,
            cuda_available = cuda_available,
            affine_yx = affine,
            velocity_shape = velocity_shape
        )
    )
}


.sa_alignment_diagnostics <- function(
    original,
    aligned,
    target,
    transformedLandmarks,
    targetLandmarks,
    maxPoints
) {
    if (!is.numeric(maxPoints) || length(maxPoints) != 1L ||
        !is.finite(maxPoints) || maxPoints < 2) {
        stop(
            "`diagnosticsMaxPoints` must be one finite number >= 2.",
            call. = FALSE
        )
    }
    bounds <- function(x) {
        list(
            x = unname(range(x[, 1L])),
            y = unname(range(x[, 2L]))
        )
    }
    output <- list(
        original_bounds = bounds(original),
        aligned_bounds = bounds(aligned),
        target_bounds = if (is.null(target)) NULL else bounds(target),
        nearest_target_before = NULL,
        nearest_target_after = NULL,
        landmark_rmse = NULL
    )
    if (!is.null(target)) {
        sample_rows <- function(n) {
            if (n <= maxPoints) seq_len(n) else unique(round(seq(1, n, length.out = maxPoints)))
        }
        target_sample <- as.matrix(target[sample_rows(nrow(target)), , drop = FALSE])
        before_sample <- as.matrix(original[sample_rows(nrow(original)), , drop = FALSE])
        after_sample <- as.matrix(aligned[sample_rows(nrow(aligned)), , drop = FALSE])
        summarize_distance <- function(query) {
            distances <- FNN::get.knnx(target_sample, query, k = 1L)$nn.dist[, 1L]
            list(
                median = unname(stats::median(distances)),
                mean = unname(mean(distances)),
                p95 = unname(stats::quantile(distances, 0.95, names = FALSE))
            )
        }
        output$nearest_target_before <- summarize_distance(before_sample)
        output$nearest_target_after <- summarize_distance(after_sample)
    }
    if (!is.null(transformedLandmarks) && !is.null(targetLandmarks)) {
        output$landmark_rmse <- sqrt(mean((transformedLandmarks - targetLandmarks)^2))
    }
    output
}
