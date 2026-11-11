#' Compute PCA on a selected experiment
#'
#' Uses [scrapper::runPca()] (IRLBA on the centred, optionally scaled, matrix
#' without materialising it) and stores embeddings in reducedDim(). Dense,
#' sparse and `DelayedArray` assays are supported.
#' The user can provide a bin/resolution size to increase the bin size and
#' reduce the dimensionality/noise of the SM dataset prior to calculating PCAs.
#'
#' @param object A SingleCellExperiment, including SpatialExperiment, containing
#'   the selected primary or alternative experiment.
#' @param npcs Positive number of PCs, capped by matrix dimensions. NULL uses
#'   the variance threshold.
#' @param varianceExplainedThreshold Cumulative variance fraction in (0, 1],
#'   used only when npcs is NULL.
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment
#'   (default = `"counts"`), read with `SummarizedExperiment::assay()`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param showVariancePlot Whether to store a variance plot in
#'   `metadata(SpaMTP)$reduction_plots[[reductionName]]` (default = FALSE).
#'   The plot is returned with the object and is never printed automatically.
#' @param binResolution Numeric value defining the resolution to use for
#'   binning m/z peaks. If set to `NULL`, no binning will be performed (default
#'   = NULL).
#' @param resolutionUnits Character string specifying the units of the
#'   `binResolution`. Either 'ppm' or 'mz' can be provided. `binResolution`
#'   must be provided for this parameter to be implemented (default = "ppm").
#' @param binMethod Character string defining the method to use for binning
#'   respective m/z peaks that fall within a bin. Options for this parameter can
#'   be one of "sum", "mean", "max" or "min". `binResolution` must be provided
#'   for this parameter to be implemented (default = "sum").
#' @param reductionName Name in `SingleCellExperiment::reducedDimNames(SpaMTP)`
#'   for the resulting embedding (default = `"pca"`). Even when the input
#'   expression is in an alternative experiment, the embedding is stored in
#'   the parent object's `reducedDim()`.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param features Optional exact feature identifiers. When supplied, all and
#'   only these rows enter PCA and `ntop` is ignored.
#' @param scale Scale features to unit variance before PCA.
#' @param ntop Number of most variable features used when `features` is NULL.
#'
#'
#' @return The input with PCA embeddings in reducedDim().
#' @export
#'
#' @examples
#' counts <- matrix(seq_len(30), nrow = 5)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#' sce <- runMetabolicPCA(sce, npcs = 2, verbose = FALSE)
#' SingleCellExperiment::reducedDimNames(sce)
#' ## For running PCA on un-adjusted peak bin sizes
#' # spamtp_obj <- runMetabolicPCA(spamtp_obj, npcs = 50)
#'
#' ## For running PCA on increased peak bin sizes
#' # spamtp_obj <- runMetabolicPCA(spamtp_obj, npcs = 50, bin_resolution = 50)
runMetabolicPCA <- function(
    object,
    npcs = 30,
    varianceExplainedThreshold = 0.9,
    assay = "main",
    assayName = "counts",
    showVariancePlot = FALSE,
    binResolution = NULL,
    resolutionUnits = "ppm",
    binMethod = "sum",
    reductionName = "pca",
    verbose = TRUE, features = NULL, scale = FALSE,
    ntop = 500
) {
    .requireExperiment(object, "SingleCellExperiment")
    analysisData <- .experimentForAssay(object, assay)
    if (!is.null(binResolution)) {
        if (!identical(analysisData, object)) {
            stop(
                "Bin the primary MSI experiment, not an alternative modality.",
                call. = FALSE
            )
        }
        analysisData <- binSpaMTP(
            analysisData,
            resolution = binResolution, units = resolutionUnits,
            method = binMethod, assay = assay, assayName = assayName
        )
        assayName <- "counts"
    }
    if (!methods::is(analysisData, "SingleCellExperiment")) {
        analysisData <- methods::as(analysisData, "SingleCellExperiment")
    }
    values <- .assayData(analysisData, layer = assayName)
    if (!is.null(features)) {
        if (!is.character(features) || anyNA(features) || anyDuplicated(features) ||
            length(features) < 2L || !all(features %in% rownames(values))) {
            stop(
                "features must contain distinct, observed feature names.",
                call. = FALSE
            )
        }
        values <- values[features, , drop = FALSE]
    }
    maximum <- min(dim(values)) - 1L
    if (maximum < 1L || !.allFiniteAssay(values)) {
        stop(
            "PCA requires finite values and at least two features and pixels.",
            call. = FALSE
        )
    }
    if (!is.null(npcs)) {
        .validateFeatureCount(npcs, "npcs")
    } else if (length(varianceExplainedThreshold) != 1L ||
        !is.finite(varianceExplainedThreshold) ||
        varianceExplainedThreshold <= 0 ||
        varianceExplainedThreshold > 1) {
        stop("varianceExplainedThreshold must be in (0, 1] when npcs is NULL.",
            call. = FALSE
        )
    }
    components <- if (is.null(npcs)) maximum else min(npcs, maximum)
    embedding <- .scrapperPca(values, components,
        scale = scale,
        ntop = if (is.null(features)) ntop else NULL
    )
    percentVar <- attr(embedding, "percentVar")
    if (is.null(npcs)) {
        retained <- which(cumsum(percentVar) / 100 >= varianceExplainedThreshold)[1L]
        if (is.na(retained)) retained <- ncol(embedding)
        loadings <- attr(embedding, "rotation")
        embedding <- embedding[, seq_len(retained), drop = FALSE]
        attr(embedding, "percentVar") <- percentVar[seq_len(retained)]
        if (!is.null(loadings)) {
            attr(embedding, "rotation") <- loadings[, seq_len(retained), drop = FALSE]
        }
    }
    if (isTRUE(showVariancePlot)) {
        varianceData <- data.frame(
            component = seq_along(percentVar), cumulative = cumsum(percentVar) / 100
        )
        S4Vectors::metadata(object)$reduction_plots[[reductionName]] <-
            ggplot2::ggplot(
                varianceData,
                ggplot2::aes(.data$component, .data$cumulative)
            ) +
            ggplot2::geom_line() +
            ggplot2::geom_point() +
            ggplot2::labs(
                x = "Principal component",
                y = "Cumulative variance explained"
            )
    }
    SingleCellExperiment::reducedDim(object, reductionName) <- embedding
    S4Vectors::metadata(object)$reduction_inputs[[reductionName]] <- list(
        assay = assay, layer = assayName, features = rownames(attr(embedding, "rotation")),
        scale = scale, method = "scrapper::runPca"
    )
    verbose_message(paste0(
        "Stored ", ncol(embedding), " principal components in reducedDim `",
        reductionName, "`."
    ), verbose = verbose)
    object
}

# PCA shared by runMetabolicPCA() and multiOmicIntegration(). With ntop, the
# most variable rows are used, matching the former scater::runPCA() default.
.scrapperPca <- function(values, ncomponents, scale = FALSE, ntop = NULL, ...) {
    if (!is.null(ntop) && ntop < nrow(values)) {
        variances <- MatrixGenerics::rowVars(DelayedArray::DelayedArray(values))
        values <- values[order(variances, decreasing = TRUE)[seq_len(ntop)], , drop = FALSE]
    }
    fit <- scrapper::runPca(values,
        number = min(ncomponents, dim(values)),
        scale = scale, ...
    )
    labels <- paste0("PC", seq_len(nrow(fit$components)))
    embedding <- t(fit$components)
    dimnames(embedding) <- list(colnames(values), labels)
    rotation <- fit$rotation
    dimnames(rotation) <- list(rownames(values), labels)
    attr(embedding, "percentVar") <- 100 * fit$variance.explained / fit$total.variance
    attr(embedding, "rotation") <- rotation
    embedding
}


#' Perform Dimensionality Reduction using Graph-Regularised PCA on Spatial Data
#'
#' Computes a graph-regularised PCA using spatial coordinates and scaled
#' expression data. A k-nearest neighbour (k-NN) graph is computed using spatial
#' locations and used to regularise the PCA decomposition via a graph Laplacian.
#' The result is stored in reducedDim(), with spatial edges in colPairs().
#'
#' Note: This method has been adapted from the GraphPCA Python package
#' (\doi{10.1186/s13059-024-03429-x}).
#'
#' @param data A SpatialExperiment or aligned Cardinal MSImagingExperiment.
#' @param nComponents Integer specifying the number of principal components to
#'   compute (default = 50).
#' @param assay Primary (`main`) or alternative experiment name.
#' @param assayName Expression assay name within the selected experiment
#'   (default = `"scaled"`), read with `SummarizedExperiment::assay()`.
#'   See [experimentAccess] for the experiment/matrix distinction.
#' @param platform Character string matching either `"Visium"` or `"ST"` to
#'   determine how the k-NN graph is constructed. If "Visium" k-nns will handle
#'   the hexagon spot arrangement, including setting `nNeighbors` = 6, else
#'   "ST" assignment will set `nNeighbors` = 4 unless a value is specifically
#'   provided (default = "Visium").
#' @param lambda Numeric value defining the regularisation parameter that
#'   controls the influence of the graph Laplacian (default = 0.5).
#' @param nNeighbors Integer value specifying the number of spatial neighbours
#'   to use. If `NULL`, will default of 6 for "Visium" data and 4 for "ST"
#'   platforms (default = NULL).
#' @param includeSelf Boolean logical value indicating whether to include
#'   self-connections in the graph (default = FALSE).
#' @param alg Character string specifying the algorithm to use for nearest
#'   neighbour search (passed to `FNN::get.knn()`) (default ="kd_tree").
#' @param fast Use partial eigendecomposition via `RSpectra::eigs_sym()`
#'   (default `TRUE`). For wide inputs requesting few components, matrix-vector
#'   products avoid forming the dense feature-by-feature covariance matrix.
#' @param graphName Name for pixel edges in `SingleCellExperiment::colPair()`
#'   and the graph matrix in `S4Vectors::metadata(data)$spatialGraphs`
#'   (default = `"SpatialKNN"`).
#' @param reductionName Name for the parent object's
#'   `SingleCellExperiment::reducedDim()` (default = `"SpatialPCA"`).
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A SpatialExperiment with a new `colPair()` graph and `reducedDim()`
#'   embedding. Expression assays are preserved.
#'
#' @details
#' The neighbour graph and Laplacian are sparse, but the selected expression
#' assay and graph-regularised solve result are materialised in memory.
#' Each dense double-precision feature-by-pixel matrix needs 8 bytes per entry:
#' 2,000 features by 100,000 pixels require about 1.6 GB per matrix, before
#' temporary copies, graph factorisation and eigensolver workspace.
#' A delayed input assay does not remove these allocations. Subset the selected
#' experiment's features or the object's pixels before calling this function
#' when the full matrices would exceed available memory. Outside the wide-input
#' operator path described by `fast`, the eigendecomposition also constructs a
#' dense feature-by-feature matrix.
#' @export
#'
#' @importFrom Matrix as.matrix colMeans colSums diag diag<- mean qr qr.resid
#' @importFrom Matrix rowMeans rowSums solve sparseMatrix summary t which
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' spe <- scaleSMData(spe)
#' spe <- runSpatialGraphPCA(spe,
#'     nComponents = 2, platform = "ST",
#'     nNeighbors = 3, verbose = FALSE
#' )
#' head(SingleCellExperiment::reducedDim(spe, "SpatialPCA"))
runSpatialGraphPCA <- function(
    data,
    nComponents = 50,
    assay = "main",
    assayName = "scaled",
    platform = "Visium",
    lambda = 0.5,
    nNeighbors = NULL,
    includeSelf = FALSE,
    alg = "kd_tree",
    fast = TRUE,
    graphName = "SpatialKNN",
    reductionName = "SpatialPCA",
    verbose = TRUE
) {
    data <- .nativeSpatialObject(data)



    if (!platform %in% c("Visium", "ST")) {
        stop("Incorrect value for plantform! platform must be assigned either 'Visium' or 'ST'... If data is not Visium 'spot'-based then set platform to 'ST'")
    }

    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        stop("lambda must be one finite, non-negative number.", call. = FALSE)
    }
    .validateFeatureCount(nComponents, "n_components")

    assayMatrix <- .assayData(data, assay, assayName)
    if (.isDelayed(assayMatrix)) {
        # The graph-regularised eigenproblem uses dense products of this matrix.
        assayMatrix <- as.matrix(assayMatrix)
    }
    if (min(dim(assayMatrix)) < 2L || !.allFiniteAssay(assayMatrix)) {
        stop("Graph PCA needs finite values and at least two features and pixels.",
            call. = FALSE
        )
    }
    nComponents <- min(nComponents, nrow(assayMatrix), ncol(assayMatrix))
    Expr <- t(assayMatrix)


    if (is.null(nNeighbors)) {
        if (platform == "Visium") {
            nNeighbors <- 6
            verbose_message(
                messageText = "No nearest neighbour value provided. Based on 'Visium' platform type setting n-neighbours = 6 ...",
                verbose = verbose
            )
        } else {
            verbose_message(
                messageText = "No nearest neighbour value provided. Based on 'ST' platform type setting n-neighbours = 4 ...",
                verbose = verbose
            )
            nNeighbors <- 4
        }
    } else {
        nNeighbors <- as.integer(nNeighbors)
        verbose_message(
            messageText = paste0(
                "Using a nearest neighbour value = ",
                nNeighbors,
                "! NOTE: the recomended values are -> 'Visium' data: n_neighbors = 6  / 'ST' data: n_neighbors = 4 ...."
            ),
            verbose = verbose
        )
    }

    location <- .nativeCoordinates(data)
    sampleIndices <- split(seq_len(nrow(location)), location$sample_id)
    graphs <- lapply(sampleIndices, function(index) {
        kNeighborsGraph(location[index, c("x", "y"), drop = FALSE],
            nNeighbors = nNeighbors, platform = platform,
            includeSelf = includeSelf, alg = alg
        )
    })
    graph <- Matrix::bdiag(graphs)
    originalOrder <- order(unlist(sampleIndices, use.names = FALSE))
    graph <- graph[originalOrder, originalOrder, drop = FALSE]
    graph <- 0.5 * (graph + t(graph))


    # Keep the Laplacian sparse: the previous conversion allocated an n-by-n
    # dense matrix and made full-tissue analysis unnecessarily expensive.
    graphL <- Matrix::Diagonal(x = Matrix::rowSums(graph)) - graph


    # Create identity matrix and add lambda * graphL
    n <- nrow(Expr)

    G <- Matrix::Diagonal(n) + (lambda * graphL)

    X <- Matrix::solve(G, as.matrix(Expr))


    rownames(graph) <- colnames(graph) <- colnames(assayMatrix)

    spatialGraphs <- S4Vectors::metadata(data)$spatialGraphs %||% list()
    spatialGraphs[[graphName]] <- graph
    S4Vectors::metadata(data)$spatialGraphs <- spatialGraphs
    edges <- summary(methods::as(graph, "generalMatrix"))
    SingleCellExperiment::colPair(data, graphName) <- S4Vectors::SelfHits(
        edges$i, edges$j,
        nnode = ncol(data), weight = edges$x
    )


    rm(G)
    rm(graph)
    rm(graphL)

    # 3. Eigendecomposition of C
    if (fast && nComponents < ncol(Expr)) {
        if (ncol(Expr) >= 256L && nComponents <= ncol(Expr) / 4) {
            # C = t(Expr) %*% X, where X = solve(G, Expr). Applying C to a vector
            # costs O(n * p), avoiding the O(n * p^2) product and O(p^2) storage.
            # The same symmetric operator and partial eigensolver are used; neither
            # observations nor features are subsampled and lambda is unchanged.
            multiply <- function(v, args) as.numeric(crossprod(Expr, X %*% v))
            eig <- RSpectra::eigs_sym(
                multiply,
                k = nComponents,
                n = ncol(Expr)
            )
        } else {
            eig <- RSpectra::eigs_sym(crossprod(Expr, X), k = nComponents)
        }
    } else {
        eig <- eigen(t(Expr) %*% X, symmetric = TRUE)
    }

    W <- eig$vectors[, seq_len(nComponents), drop = FALSE]
    Z <- X %*% W

    rownames(Z) <- colnames(assayMatrix) # cells
    rownames(W) <- rownames(assayMatrix) # genes (or features)
    colnames(Z) <- paste0("PC_", seq_len(ncol(Z)))
    colnames(W) <- paste0("PC_", seq_len(ncol(W)))

    # Create the PCA reduction
    SingleCellExperiment::reducedDim(data, reductionName) <- as.matrix(Z)
    loadings <- S4Vectors::metadata(data)$reductionLoadings %||% list()
    loadings[[reductionName]] <- W
    S4Vectors::metadata(data)$reductionLoadings <- loadings
    S4Vectors::metadata(data)$reduction_inputs[[reductionName]] <- list(
        assay = assay, layer = assayName, features = rownames(assayMatrix),
        observations = colnames(assayMatrix), lambda = lambda, neighbors = nNeighbors,
        graph = graphName, sample_separated = TRUE, method = "runSpatialGraphPCA"
    )


    return(data)
}


#' Construct a k-Nearest Neighbour Graph from Spatial Coordinates.
#'
#' Builds a sparse adjacency matrix representing the k-nearest neighbour
#' relationships between spots or cells, using Euclidean distance between
#' spatial coordinates.
#'
#' @param location A data frame or matrix with columns `x` and `y` representing
#'   spatial coordinates, and rownames match barcode/cell names.
#' @param nNeighbors Integer value defining the number of neighbours to
#'   consider.
#' @param platform Character string being either `"Visium"` or `"ST"`,
#'   specifying the platform of the dataset being analysed. Determines how
#'   distance thresholds and adjacency are calculated whereby "Visium" data is
#'   handled slightly differently to compensate for the hexagon shape of the
#'   spots.
#' @param includeSelf Boolean logical value stating whether to include
#'   self-connections (default = FALSE).
#' @param alg Character string specifying the algorithm to use for nearest
#'   neighbour search. Passed to `FNN::get.knn()` (default = "kd_tree").
#'
#' @return A sparse adjacency matrix (`dgCMatrix`) representing the k-NN graph.
#' @export
#'
#' @examples
#' coordinates <- cbind(x = rep(1:3, 2), y = rep(1:2, each = 3))
#' kNeighborsGraph(coordinates, nNeighbors = 2, platform = "ST")
kNeighborsGraph <- function(
    location,
    nNeighbors,
    platform,
    includeSelf = FALSE,
    alg = "kd_tree"
) {
    .validateFeatureCount(nNeighbors, "n_neighbors")
    if (any(!is.finite(as.matrix(location)))) {
        stop("Spatial coordinates must be finite.", call. = FALSE)
    }
    if (nrow(location) <= 1L) {
        return(Matrix::Diagonal(nrow(location), x = as.numeric(includeSelf)))
    }
    nNeighbors <- min(nNeighbors, nrow(location) - 1L)
    # Get the k-nearest neighbors using kd_tree (most similar to scikit-learn's default)
    knn_result <- FNN::get.knn(location, k = nNeighbors, algorithm = alg)

    if (platform == "Visium") {
        # Step 1: Flatten the distances and indices
        N <- nrow(knn_result$nn.index)

        # Reshape distances and column indices
        dists <- as.vector(t(knn_result$nn.dist))
        col_indices <- as.vector(t(knn_result$nn.index))

        # Create row indices (equivalent to np.repeat())
        row_indices <- rep(seq_len(N), each = nNeighbors)

        # Apply neighbor correction if needed

        dist_cutoff <- median(dists) * 1.3 # Small amount of sway
        mask <- dists < dist_cutoff

        # Filter using the mask
        row_indices <- row_indices[mask]
        col_indices <- col_indices[mask]
        dists <- dists[mask]

        adjacency <- sparseMatrix(
            i = row_indices,
            j = col_indices,
            x = rep(1, length(row_indices)),
            dims = c(N, N),
            repr = "C" # Ensures a "dgCMatrix" (like CSR in Python)
        )
    } else {
        # Create empty adjacency matrix
        n <- nrow(location)
        adjacency <- Matrix::sparseMatrix(
            i = rep(seq_len(n), each = nNeighbors),
            j = as.vector(t(knn_result$nn.index)), x = 1, dims = c(n, n)
        )
    }

    # Include self-loops if requested (in this case, no)
    if (includeSelf) {
        diag(adjacency) <- 1
    }


    return(adjacency)
}


#' Perform K-means clustering on a specified reduction
#'
#' Runs K-means on an embedding read through
#' `SingleCellExperiment::reducedDim()`
#' and writes pixel cluster labels to `SummarizedExperiment::colData()`.
#'
#' @param data A SingleCellExperiment, including SpatialExperiment, containing
#'   an embedding such as the result of `runMetabolicPCA()`.
#' @param reduction Name in `SingleCellExperiment::reducedDimNames(data)`
#'   (default = `"pca"`).
#' @param clusterName Column name in `SummarizedExperiment::colData(data)`
#'   for cluster labels (default = `"spatial_clusters"`).
#' @param clusters Integer defining the number of clusters to form (default =
#'   8).
#' @param iterMax Integer defining the maximum number of iterations allowed
#'   (default = 10).
#' @param nstart Integer stating the number of random sets to choose (default =
#'   1).
#' @param algorithm Character string defining the K-means algorithm to use. One
#'   of `"Hartigan-Wong"`, `"Lloyd"`, `"Forgy"`, or `"MacQueen"` (default =
#'   "Hartigan-Wong").
#' @param trace Logical boolean indicating whether to produce tracing
#'   information on the progress of the algorithm (default = FALSE).
#' @param seed Integer of the random seed to use for reproducibility (default =
#'   888).
#'
#' @return The input container with cluster assignments in
#'   `SummarizedExperiment::colData(data)[[clusterName]]`.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' spe <- getKmeanClusters(spe, reduction = "pca", clusters = 2)
#' table(spe$spatial_clusters)
getKmeanClusters <- function(
    data,
    reduction = "pca",
    clusterName = "spatial_clusters",
    clusters = 8,
    iterMax = 10,
    nstart = 1,
    algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
    trace = FALSE,
    seed = 888
) {
    .requireExperiment(data, "SingleCellExperiment")


    if (!reduction %in% SingleCellExperiment::reducedDimNames(data)) {
        stop(
            "Reduction `",
            reduction,
            "` is not present in the object.",
            call. = FALSE
        )
    }
    embedding <- SingleCellExperiment::reducedDim(data, reduction)

    labels <- .kmeanLabels(
        embedding,
        clusters,
        iterMax,
        nstart,
        algorithm,
        trace,
        seed
    )

    metadata <- .cellMetadata(data)
    metadata[[clusterName]] <- labels
    data <- .setCellMetadata(data, metadata)
    return(data)
}

.kmeanLabels <- function(
    embedding, clusters, iterMax = 100, nstart = 20,
    algorithm = "Hartigan-Wong", trace = FALSE, seed = 888
) {
    withr::with_seed(
        seed,
        stats::kmeans(
            embedding,
            centers = clusters,
            iter.max = iterMax,
            nstart = nstart,
            algorithm = algorithm,
            trace = trace
        )$cluster
    )
}
