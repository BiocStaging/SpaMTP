# Download only images with URLs; missing isotopes remain explicit NULL entries.
.metaspace_image_pixels <- function(imgInfo, relative = TRUE, hotspotClipping = FALSE) {
    if (is.null(imgInfo$url) || !nzchar(imgInfo$url)) return(NULL)
    # ---- download PNG
    response <- httr::GET(imgInfo$url)
    httr::stop_for_status(response)
    img_raw <- httr::content(response, "raw")
    img <- png::readPNG(img_raw)
    nchan <- if (length(dim(img)) == 2L) 1L else dim(img)[3]

    # ---- intensity channel
    # R PNG reads array as [row, col, channel]
    if (length(dim(img)) == 2L) {
        intensity_data <- img
    } else if (nchan %in% c(1, 2, 3, 4)) {
        intensity_data <- matrix(img[, , 1], nrow = dim(img)[1], ncol = dim(img)[2])
    } else {
        stop("Unsupported PNG: ", nchan, " channels")
    }

    # ---- scale to real intensity
    lo <- imgInfo$minIntensity %||% 0
    hi <- imgInfo$maxIntensity %||% 1

    # Perform scaling (0-1 normalized to min/max intensity)
    intensities <- lo + intensity_data * (hi - lo)

    # ---- alpha mask (Apply mask to set non-data points to 0)
    if (nchan %in% c(2, 4)) {
        mask <- img[, , nchan]
        intensities[mask == 0] <- 0
    }

    # ---- relative / absolute (Normalization)
    if (relative && hi > 0) intensities <- intensities / hi
    # If relative=FALSE, intensities are absolute (raw counts).

    # ---- hotspot clipping (Always apply to the final intensity values)
    if (hotspotClipping && hi > 0) { # Only clip if there's actual signal
        # Match the Python client's 99th percentile above the visible floor.
        positive_intensities <- intensities[intensities > max(intensities) / 256]
        if (length(positive_intensities) > 0) {
            limit <- stats::quantile(positive_intensities, 0.99, type = 7)
            intensities[intensities > limit] <- limit
        }
    }

    # ---- RETURN FULL INTENSITY MATRIX (2D array)
    return(intensities)
}

#### SpaMTP METASPACE R CLIENT #################################################################################################################################################################################


#' @title Metaspace R Client Constructor
#'
#' @description Creates a client object with methods to interact with the
#'   METASPACE GraphQL API.
#'
#' @param host Character string specifying the METASPACE host URL (default = "https://metaspace2020.org").
#' @param apiKey Optional. Your METASPACE API key for accessing private
#'   datasets. If accessing a public dataset this can be left as `NULL` (default
#'   = NULL).
#' @param BPPARAM BiocParallel parameters for downloading annotation images.
#'   Defaults to serial execution without changing the session's parallel plan.
#'
#' @return A list containing functions for querying METASPACE data. Annotation
#'   results preserve `isotopeImages` as a list column. Missing isotope URLs
#'   yield NULL image entries. Hotspot clipping follows the Python client:
#'   the 99th percentile of pixels brighter than the image maximum / 256.
#'
#' @importFrom httr add_headers POST GET status_code content
#' @importFrom jsonlite toJSON
#' @importFrom rlang %||%
#'
#' @export
#'
#' @examples
#' client <- metaspaceClient(BPPARAM = BiocParallel::SerialParam())
#' names(client)
#' # Querying the remote service requires an internet connection.
#' if (interactive()) client$get_info("2020-12-07_03h16m14s")
metaspaceClient <- function(
    host = "https://metaspace2020.org", apiKey = NULL,
    BPPARAM = BiocParallel::SerialParam()
) {
    graphql_url <- paste0(host, "/graphql")


    ##  GraphQL executor
    query_graphql <- function(query, variables = list()) {
        hdr <- list(`Content-Type` = "application/json", `Source` = "api")
        if (!is.null(apiKey)) hdr$Authorization <- paste0("Api-Key ", apiKey)
        headers <- do.call(add_headers, hdr)

        resp <- POST(graphql_url,
            body = toJSON(
                list(query = query, variables = variables),
                auto_unbox = TRUE
            ),
            headers, encode = "json"
        )
        if (status_code(resp) == 401) {
            stop("Authentication required - get your key at https://metaspace2020.org/user/me")
        }
        if (status_code(resp) != 200) stop("HTTP ", status_code(resp))
        res <- content(resp, as = "parsed")
        if (!is.null(res$errors)) {
            msg <- res$errors[[1]]$message
            stop("GraphQL request failed: ", msg)
        }
        res$data
    }


    ##  1. annotation table (with pagination)
    get_dataset_results <- function(
        datasetId, fdr = 0.1,
                                            database = c("HMDB", "v4"),
                                            includeChemMods = FALSE,
                                            includeNeutralLosses = FALSE,
                                            pageSize = 100
    ) {
        # ---- get DB ID
        db_res <- query_graphql("query { allMolecularDBs { id name version } }")
        db_id <- NULL
        for (db in db_res$allMolecularDBs) {
            if (db$name == database[1] && db$version == database[2]) {
                db_id <- db$id
                break
            }
        }
        if (is.null(db_id)) stop(
            "Database not found: ",
            paste(database, collapse = " ")
        )

        # ---- construct filter
        filter <- list(fdrLevel = fdr, databaseId = db_id)
        if (!includeChemMods) filter$hasChemMod <- FALSE
        if (!includeNeutralLosses) filter$hasNeutralLoss <- FALSE

        # ---- GraphQL query with pagination
        q <- "
    query($filter: AnnotationFilter, $dFilter: DatasetFilter, $limit: Int, $offset: Int) {
      allAnnotations(filter: $filter, datasetFilter: $dFilter, limit: $limit, offset: $offset) {
        sumFormula
        adduct
        neutralLoss
        chemMod
        mz
        msmScore
        rhoSpatial
        rhoSpectral
        rhoChaos
        fdrLevel
        offSample
        possibleCompounds { name }
        isotopeImages { mz url minIntensity maxIntensity totalIntensity }
      }
    }
  "

        ann_list <- list()
        offset <- 0
        repeat {
            vars <- list(
                filter = filter, dFilter = list(ids = datasetId),
                limit = pageSize, offset = offset
            )
            res <- query_graphql(q, vars)

            anns <- res$allAnnotations
            if (length(anns) == 0) break # no more annotations

            # ---- convert annotations to data.frame
            page_df <- lapply(anns, function(a) {
                names_str <- paste(
                    vapply(a$possibleCompounds, function(compound) compound$name %||% "", character(1)),
                    collapse = "; "
                )
                int <- if (length(a$isotopeImages)) a$isotopeImages[[1]]$maxIntensity else NA_real_
                row <- data.frame(
                    formula = a$sumFormula, adduct = a$adduct,
                    chemMod = ifelse(is.null(a$chemMod), "", a$chemMod),
                    neutralLoss = ifelse(is.null(a$neutralLoss), "", a$neutralLoss),
                    mz = a$mz, msm = a$msmScore %||% NA_real_, fdr = a$fdrLevel,
                    rhoSpatial = a$rhoSpatial %||% NA_real_, rhoSpectral = a$rhoSpectral %||% NA_real_,
                    moc = a$rhoChaos %||% NA_real_, offSample = a$offSample %||% NA,
                    intensity = int %||% NA_real_, moleculeNames = names_str,
                    stringsAsFactors = FALSE
                )
                row$isotopeImages <- list(a$isotopeImages)
                row
            })
            ann_list <- c(ann_list, page_df)

            offset <- offset + length(anns)
        }

        if (length(ann_list) == 0) {
            warning("No annotations found")
            return(data.frame())
        }

        do.call(rbind, ann_list)
    }


    ##  2. dataset metadata
    get_dataset_info <- function(datasetId) {
        q <- "
      query($id:String!){dataset(id:$id){
        id name status uploadDT polarity ionisationSource
        analyzer{type resolvingPower(mz:400)}
        organism organismPart condition maldiMatrix isPublic
        adducts submitter{id name} group{id name}
        databases{id name version}
      }}"
        res <- query_graphql(q, list(id = datasetId))
        if (is.null(res$dataset)) stop("Dataset not found")
        res$dataset
    }


    ##  3. PNG dimension helper
    .png_dim <- function(url) {
        if (status_code(HEAD(url)) != 200) stop("Cannot reach image")
        raw <- GET(url, config = list(range = "bytes=0-23"))$content
        w <- readBin(raw[17:20], "integer", size = 4, endian = "big")
        h <- readBin(raw[21:24], "integer", size = 4, endian = "big")
        c(width = w, height = h)
    }


    ##  4. robust pixel extractor (MODIFIED TO RETURN FULL MATRIX)
    get_ion_image_pixels <- function(
        datasetId, sumFormula, adduct,
                                             fdr = 0.1, database = c("HMDB", "v4"),
                                             chemMod = NULL, neutralLoss = NULL,
                                             isotopeIdx = 1, relative = TRUE,
                                             hotspotClipping = FALSE
    ) {
        # ---- DB id
        db_res <- query_graphql("query{allMolecularDBs{id name version}}")
        db_id <- NULL
        for (db in db_res$allMolecularDBs) {
            if (db$name == database[1] && db$version == database[2]) {
                db_id <- db$id
                break
            }
        }
        if (is.null(db_id)) stop("Database not found")

        # ---- filter
        filter <- list(
            sumFormula = sumFormula, adduct = adduct,
            databaseId = db_id, fdrLevel = fdr
        )
        if (!is.null(chemMod) && chemMod != "") filter$chemMod <- chemMod
        if (!is.null(neutralLoss) && neutralLoss != "") filter$neutralLoss <- neutralLoss

        q <- "
      query($filter:AnnotationFilter!,$dFilter:DatasetFilter){
        allAnnotations(filter:$filter,datasetFilter:$dFilter){
          isotopeImages{url mz maxIntensity minIntensity totalIntensity}
        }}"
        vars <- list(filter = filter, dFilter = list(ids = datasetId))
        res <- query_graphql(q, vars)
        if (!length(res$allAnnotations)) stop("Annotation not found")
        available <- res$allAnnotations[[1]]$isotopeImages
        if (length(isotopeIdx) != 1L || !is.finite(isotopeIdx) ||
            isotopeIdx < 1 || isotopeIdx != floor(isotopeIdx)) {
            stop("isotopeIdx must be a positive integer.", call. = FALSE)
        }
        if (isotopeIdx > length(available) || is.null(available[[isotopeIdx]]$url)) {
            return(NULL)
        }
        img_info <- available[[isotopeIdx]]

        .metaspace_image_pixels(img_info, relative, hotspotClipping)
    }


    ##  5. Python-style all_annotation_images
    all_annotation_images <- function(
        datasetId,
                                              fdr = 0.1,
                                              database = c("HMDB", "v4"),
                                              onlyFirstIsotope = TRUE,
                                              relative = TRUE,
                                              hotspotClipping = FALSE
    ) {
        ann_df <- get_dataset_results(
            datasetId,
            fdr = fdr,
            database = database
        )
        if (nrow(ann_df) == 0) stop("No annotations found at FDR ", fdr)

        get_single_image <- function(row) {
            isotope_idx <- if (onlyFirstIsotope) 1 else seq_along(row$isotopeImages[[1]])

            lapply(isotope_idx, function(idx) {
                images <- row$isotopeImages[[1]]
                img_info <- if (idx <= length(images)) images[[idx]] else NULL
                img_df <- .metaspace_image_pixels(img_info, relative, hotspotClipping)
                if (is.null(img_df) || !nrow(img_df)) {
                    return(NULL)
                }
                img_df
            })
        }

        all_images <- BiocParallel::bplapply(seq_len(nrow(ann_df)),
            function(i) get_single_image(ann_df[i, , drop = FALSE]),
            BPPARAM = BPPARAM
        )
        names(all_images) <- ann_df$mz
        all_images
    }

    ##  Return public API (unchanged)
    list(
        query                 = query_graphql,
        get_results           = get_dataset_results,
        get_info              = get_dataset_info,
        get_ion_image_pixels  = get_ion_image_pixels,
        all_annotation_images = all_annotation_images,
        .png_dim              = .png_dim
    )
}


#' @title Core METASPACE Data Retrieval
#'
#' @description The main user-facing function to download annotations and ion
#'   images for a METASPACE dataset.
#'
#' @param datasetId Character string of the METASPACE dataset ID.
#' @param fdr Numeric, the FDR level threshold (default = 0.1).
#' @param database Character vector, specifying the database name and version
#'   (default = c("HMDB", "v4")).
#' @param includeImages Logical, if TRUE, downloads all ion image matrices
#'   (default = TRUE).
#' @param apiKey Optional character string, your METASPACE API key (default =
#'   NULL).
#' @param isotopeIdx Integer, index of the isotope image to download for all
#'   annotations (1 for main peak) (default = 1).
#' @param relative Logical, if TRUE, normalizes pixel intensities (0-1 range)
#'   (default = TRUE).
#' @param verbose Logical; show download messages and a progress bar (default =
#'   TRUE).
#'
#' @return A list with `annotations` (data.frame) and `images` (list of 2D
#'   matrices, with NULL for unavailable images). Feature names are unique
#'   even when annotations have the same rounded m/z; the numeric `mz` column
#'   retains full precision.
#'
#' @importFrom httr HEAD status_code
#' @inheritParams metaspaceClient
#' @importFrom dplyr select
#' @export
#'
#' @examples
#' # Requires the METASPACE service and an internet connection.
#' if (interactive()) {
#'     result <- getMetaspace("2020-12-07_03h16m14s",
#'         includeImages = FALSE,
#'         verbose = FALSE
#'     )
#'     head(result$annotations)
#' }
getMetaspace <- function(
    datasetId,
    fdr = 0.1,
    database = c("HMDB", "v4"),
    includeImages = TRUE,
    apiKey = NULL,
    isotopeIdx = 1,
    relative = TRUE,
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam(progressbar = verbose)
) {
    sm <- metaspaceClient(apiKey = apiKey)
    if (!is.numeric(isotopeIdx) || length(isotopeIdx) != 1L ||
        !is.finite(isotopeIdx) || isotopeIdx < 1 || isotopeIdx != floor(isotopeIdx)) {
        stop("isotopeIdx must be a positive integer.", call. = FALSE)
    }

    ann_df <- sm$get_results(datasetId, fdr = fdr, database = database)
    if (nrow(ann_df) == 0) stop("No annotations at FDR ", fdr)
    ann_df <- ann_df[order(ann_df$mz), ]
    mz_names <- make.unique(sprintf("%.5f", ann_df$mz))
    rownames(ann_df) <- mz_names

    images <- NULL
    if (includeImages) {
        verbose_message(
            paste0(
                "Downloading ",
                nrow(ann_df),
                " ion images (full matrices) ..."
            ),
            verbose = verbose
        )
        mats <- BiocParallel::bplapply(seq_len(nrow(ann_df)), function(i) {
            row <- ann_df[i, , drop = FALSE]
            available <- row$isotopeImages[[1]]
            img_info <- if (isotopeIdx <= length(available)) available[[isotopeIdx]] else NULL
            .metaspace_image_pixels(img_info, relative = relative)
        }, BPPARAM = BPPARAM)
        names(mats) <- mz_names
        images <- mats
    }

    list(annotations = ann_df, images = images)
}


#' @title Convert METASPACE Images to Feature Matrix
#'
#' @description Transforms the list of 2D ion image matrices from \code{getMetaspace} into a single matrix suitable for spatial analysis packages.
#'
#' @param metaspaceData A list object returned by \code{getMetaspace}.
#' @param transform Logical. If TRUE, transposes each image matrix before
#'   flattening (e.g., to align coordinate systems) (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A numeric matrix where rows are molecular features (m/z) and columns
#'   are spatial pixels (named 'x_y'). Missing images are omitted with a warning;
#'   if none are available, the function stops.
#'
#' @export
#'
#' @examples
#' images <- list(`100.00000` = matrix(1:6, 2), `200.00000` = matrix(6:1, 2))
#' metaspaceToFeatureMatrix(list(images = images), verbose = FALSE)
metaspaceToFeatureMatrix <- function(
    metaspaceData,
    transform = FALSE,
    verbose = TRUE
) {
    if (is.null(metaspaceData$images) || length(metaspaceData$images) == 0) {
        stop("Input 'metaspaceData' does not contain image matrices. Ensure includeImages=TRUE was used.")
    }
    missing_images <- vapply(metaspaceData$images, is.null, logical(1))
    if (all(missing_images)) {
        stop("No downloadable ion images are available for this selection.", call. = FALSE)
    }
    if (any(missing_images)) {
        warning("Omitting ", sum(missing_images), " annotation(s) without an ion image from the feature matrix.",
            call. = FALSE)
        metaspaceData$images <- metaspaceData$images[!missing_images]
    }
    if (transform) {
        metaspaceData$images <- lapply(metaspaceData$images, t)
    }

    # 1. Extract first image to get dimensions and create column names
    first_image <- metaspaceData$images[[1]]


    if (!is.matrix(first_image)) {
        stop("Image data must be a matrix (full pixel grid). Ensure get_ion_image_pixels returns a matrix.")
    }

    # Get dimensions (rows = y, cols = x)
    h <- nrow(first_image)
    w <- ncol(first_image)

    # Create column names (pixel coordinates 'x_y')
    # We iterate x from 1 to w, and y from 1 to h
    x_coords <- rep(seq_len(w), each = h)
    y_coords <- rep(seq_len(h), times = w)
    pixel_names <- paste0(x_coords, "_", y_coords)

    # 2. Flatten all image matrices and combine into the final matrix
    # The final matrix will be (Features x Pixels)

    verbose_message(paste0(
        "Combining ", length(metaspaceData$images),
        " annotations into a feature matrix..."
    ), verbose)

    # Check for consistency
    consistent <- vapply(metaspaceData$images, function(mat) {
        is.matrix(mat) && nrow(mat) == h && ncol(mat) == w
    }, logical(1))
    if (!all(consistent)) {
        stop(
            "Image matrix ",
            which(!consistent)[1],
            " has inconsistent dimensions."
        )
    }

    flat_images_list <- lapply(metaspaceData$images, as.vector)

    # Combine list of vectors into a matrix (Features x Pixels)
    feature_matrix <- do.call(rbind, flat_images_list)

    # Assign column names (pixel coordinates)
    colnames(feature_matrix) <- pixel_names

    verbose_message(
        messageText = paste0(
            "Conversion complete. Matrix dimensions: ",
            paste(dim(feature_matrix), collapse = " rows x "),
            " columns"
        ),
        verbose = verbose
    )

    return(feature_matrix)
}


#' @title Load METASPACE Data as a SpatialExperiment
#'
#' @description Downloads METASPACE data, converts it to a feature matrix, and
#' packages it in a Bioconductor `SpatialExperiment`. Use explicit conversion
#' for other containers.
#' Coordinates are stored in `spatialCoords()` without redundant copies in
#' `colData()`; dataset and acquisition metadata remain in `colData()`.
#'
#' @param datasetId Character string of the METASPACE dataset ID.
#' @param fdr Numeric, the FDR level threshold (default = 0.1).
#' @param database Character vector, specifying the database name and version
#'   (default = c("HMDB", "v4")).
#' @param apiKey Optional character string, your METASPACE API key (default =
#'   NULL).
#' @param relative Logical, if TRUE, normalizes pixel intensities (0-1 range)
#'   (default = FALSE).
#' @param transform Boolean specifiying whether to transform/flip the image
#'   (default = FALSE).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param returnType Output container; only `"SpatialExperiment"` is supported.
#' @inheritParams getMetaspace
#'
#' @return A `SpatialExperiment`.
#'
#' @importFrom tidyr separate
#' @importFrom png readPNG
#' @importFrom dplyr select
#'
#' @export
#'
#' @examples
#' # Requires the METASPACE service and downloads image data.
#' if (interactive()) {
#'     spe <- loadMetaspace("2020-12-07_03h16m14s", verbose = FALSE)
#'     plotSpatialFeature(spe, rownames(spe)[1])
#' }
loadMetaspace <- function(
    datasetId,
    fdr = 0.1,
    database = c("HMDB", "v4"),
    apiKey = NULL,
    relative = FALSE,
    transform = FALSE,
    verbose = TRUE,
    returnType = "SpatialExperiment",
    BPPARAM = BiocParallel::SerialParam(progressbar = verbose)
) {
    returnType <- match.arg(returnType)

    verbose_message(
        messageText = paste0(
            "Downloading dataset [",
            datasetId,
            "] from METASPACE... "
        ),
        verbose = verbose
    )

    data <- getMetaspace(
        datasetId = datasetId, fdr = fdr, database = database,
        includeImages = TRUE, apiKey = apiKey, relative = relative,
        verbose = verbose, BPPARAM = BPPARAM
    )

    verbose_message(
        messageText = "Gathering metabolite intensity values per pixel ...",
        verbose = verbose
    )

    sparse_matrix <- metaspaceToFeatureMatrix(
        data,
        transform = transform,
        verbose = verbose
    )
    features <- rownames(sparse_matrix)
    pixels <- colnames(sparse_matrix)

    pixel_data <- data.frame(pixel = pixels) %>% separate(
        col = "pixel",
        into = c("x", "y"),
        sep = "_",
        convert = TRUE
    )

    featureNames <- make.unique(paste0("mz-", features))
    rownames(sparse_matrix) <- featureNames
    pixel_metadata <- S4Vectors::DataFrame(
        metaspace_dataset = rep(datasetId, length(pixels)),
        fdr_threshold = rep(fdr, length(pixels)),
        database = rep(paste(database, collapse = "-"), length(pixels)),
        row.names = pixels
    )

    metaspace_metadata_names <- c(
        "formula", "adduct", "chemMod", "neutralLoss", "mz", "msm",
        "fdr", "rhoSpatial", "rhoSpectral", "moc", "offSample", "intensity", "moleculeNames"
    )

    availableMetadata <- intersect(
        metaspace_metadata_names,
        colnames(data$annotations)
    )
    annotation_rows <- match(features, names(data$images))
    if (anyNA(annotation_rows) || any(annotation_rows > nrow(data$annotations))) {
        stop("Ion image names cannot be aligned to annotation rows.", call. = FALSE)
    }
    featureMetadata <- data$annotations[annotation_rows, availableMetadata, drop = FALSE]
    featureMetadata$raw_mz <- featureMetadata$mz
    featureMetadata$mz <- featureMetadata$raw_mz
    featureMetadata$mz_names <- featureNames
    rownames(featureMetadata) <- featureNames
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = Matrix::Matrix(sparse_matrix, sparse = TRUE)),
        rowData = S4Vectors::DataFrame(featureMetadata),
        colData = S4Vectors::DataFrame(pixel_metadata),
        spatialCoords = as.matrix(pixel_data[, c("x", "y"), drop = FALSE]),
        sample_id = datasetId
    )
    S4Vectors::metadata(object)$metaspace <- list(
        dataset_id = datasetId,
        fdr = fdr,
        database = database,
        relative = relative
    )
    object
}
