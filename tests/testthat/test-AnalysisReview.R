test_that("unequal pixel populations do not create technical fold changes", {
    withr::local_seed(90)
    x <- matrix(stats::rpois(40 * 3300, 100), 40,
        dimnames = list(paste0("f", 1:40), paste0("p", 1:3300)))
    object <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = x),
        colData = S4Vectors::DataFrame(region = factor(rep(c("large", "small"), c(3000, 300)))))
    expect_warning(null <- findAllDEMs(object, "region", method = "technical", verbose = FALSE), "technical subdivisions")
    expect_true(all(null$DEMs$FDR > 0.05))
    expect_lt(max(abs(null$DEMs$logFC)), 0.05)
    SummarizedExperiment::assay(object)[1:4, 1:3000] <- x[1:4, 1:3000] * 2
    expect_warning(signal <- findAllDEMs(object, "region", method = "technical", verbose = FALSE), "technical subdivisions")
    hit <- signal$DEMs[signal$DEMs$cluster == "large" & signal$DEMs$gene %in% paste0("f", 1:4), ]
    expect_equal(nrow(hit), 4L)
    expect_true(all(hit$FDR < 0.05 & hit$logFC > 0.9 & hit$logFC < 1.1))
    default <- findAllDEMs(object, "region", verbose = FALSE)
    expect_true(all(c("mean_auc", "mean_cohen") %in% names(default$DEMs)))
    expect_false(any(c("FDR", "p_val_adj", "P.Value") %in% names(default$DEMs)))
})

test_that("numeric spectral tolerances have explicit mass units", {
    x <- matrix(1:6, 3, dimnames = list(NULL, c("a", "b")))
    absolute <- spectralBinning(x, c(100, 200), c(100, 100.01, 200), tolerance = .02)
    expect_equal(unname(absolute), rbind(c(3, 9), c(3, 6)))
    expect_identical(colnames(absolute), colnames(x))
    expect_equal(absolute, spectralBinning(x, c(100, 200), c(100, 100.01, 200), tolerance = c(relative = 2e-4)))
    expect_error(spectralBinning(x, 100, 1:2, tolerance = .02), "per matrix row")
    expect_error(spectralBinning(x, 100, 1:3, tolerance = NA_real_), "positive finite")
})

test_that("unnamed observations work through PCA, clustering and ion plots", {
    withr::local_seed(17)
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = matrix(stats::runif(30), 3, dimnames = list(c("a", "b", "c"), NULL))),
        rowData = S4Vectors::DataFrame(mz = c(100, 200, 300)),
        spatialCoords = cbind(x = 1:10, y = rep(1:2, 5)))
    before <- object
    object <- runMetabolicPCA(object, npcs = 2, verbose = FALSE)
    clustered <- getKmeanClusters(object, clusters = 2)
    expect_equal(length(clustered$spatial_clusters), 10L)
    expect_null(colnames(clustered))
    plot <- imageMZPlot(before, 100)
    expect_s3_class(plot, "ggplot")
    expect_equal(plot$data$score, as.numeric(SummarizedExperiment::assay(before)[1, ]))
    expect_null(colnames(before))
})

test_that("finite validation does not allocate a dense copy of a sparse assay", {
    x <- Matrix::sparseMatrix(i = seq_len(5000), j = rep(1:1000, 5), x = 1, dims = c(5000, 1000))
    expect_true(.allFiniteAssay(x))
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
        bad <- x
        bad[1, 1] <- value
        expect_false(.allFiniteAssay(bad))
    }
    if (capabilities("profmem")) {
        path <- withr::local_tempfile()
        utils::Rprofmem(path)
        .allFiniteAssay(x)
        utils::Rprofmem(NULL)
        bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(path))))
        expect_lt(max(c(0, bytes), na.rm = TRUE), 8e6)
    }
})

test_that("regional pathway selection works without a source ID column", {
    pathways <- data.frame(pathwayRampId = c("P1", "P2"), pathwayName = c("one", "two"),
        Cluster_id = "A", NES = c(2, -1), size = 3, padj = .01)
    db <- list(pathway = pathways[, 1:2], analytehaspathway = data.frame(pathwayRampId = c("P1", "P2"), rampId = c("C1", "C2")))
    plot <- plotRegionalPathways(pathways, selected_pathways = "P1", database = db)
    expect_identical(unique(plot$data$pathwayRampId), "P1")
    expect_error(plotRegionalPathways(pathways, selected_pathways = "absent", database = db), "No selected pathways")
})

test_that("pathway thumbnails handle small, constant and irregular panels", {
    coords <- data.frame(x = c(0, 2, 0), y = c(0, 0, 2))
    for (n in 1:2) {
        image <- .pathwayFeatureRaster(matrix(1, 3, n), coords)
        expect_s3_class(image, "raster")
        expect_equal(dim(image), c(2L, 2L))
        expect_equal(sum(grDevices::col2rgb(as.vector(image), alpha = TRUE)[4, ] == 0), 1L)
    }
    expect_s3_class(.pathwayFeatureRaster(matrix(1), coords[1, ]), "raster")
})

test_that("pathway dendrograms cluster shared membership first", {
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = matrix(1:6, 1, dimnames = list("mz-100.1", paste0("p", 1:6)))),
        spatialCoords = cbind(x = 1:6, y = 0))
    pathways <- data.frame(pathway_name = c("one", "two", "three"), pathway_id = c("S1", "S2", "S3"),
        p_val = .01, analytes_in_pathways = 1L, total_in_pathways = 3L, adduct_info = "mz-100.1")
    db <- list(pathway = data.frame(pathwayRampId = c("P1", "P2", "P3"), sourceId = c("S1", "S2", "S3")),
        analytehaspathway = data.frame(pathwayRampId = c("P1", "P1", "P2", "P2", "P3"), rampId = c("A", "B", "A", "B", "C")))
    distances <- NULL
    testthat::local_mocked_bindings(hclust = function(d, ...) {
        distances <<- as.matrix(d)
        stats::hclust(d, ...)
    })
    expect_s3_class(visualisePathways(object, pathways, assay = "main", min_n = 1,
        database = db, verbose = FALSE), "ggplot")
    expect_equal(unname(distances), matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3))
})

test_that("normalization has the same log2 units in both interfaces", {
    data(spamtp_example)
    object <- spamtp_example$experiment
    normalized <- normalizeSMData(object, "LogNormalize", verbose = FALSE)
    expected <- .wf_normalize(SummarizedExperiment::assay(object, "counts"), list(normalization = "library_log", scale_factor = 1e4))
    expect_equal(SummarizedExperiment::assay(normalized, "logcounts"), expected)
    expect_identical(S4Vectors::metadata(normalized)$spamtp_normalization$log_base, 2)
    summary <- summarizeFeatures(list(first = object, second = object), features = rownames(object)[1], groupBy = "region")
    expect_equal(summary$mean[1:2], summary$mean[3:4])
    expect_equal(sum(summary$n), ncol(object) * 2)
    expect_error(pathwayMapping(object), "No pathway mapping")
    expect_error(workflowMarkers(list()), "No regional markers")
})

test_that("named differential tables are independent of list ordering", {
    data(spamtp_example)
    object <- spamtp_example$experiment
    db <- spamtp_example$database
    genes <- data.frame(gene = rownames(SingleCellExperiment::altExp(object, "transcriptome")),
        cluster = "left", avg_log2FC = seq_len(8), p_val_adj = .01)
    metabolites <- data.frame(gene = rownames(object), cluster = "left",
        avg_log2FC = seq_len(nrow(object)), p_val_adj = .01)
    seen <- list()
    testthat::local_mocked_bindings(fgsea = function(pathways, stats, ...) {
        seen[[length(seen) + 1L]] <<- stats
        data.table::data.table(pathway = character(), pval = numeric(), padj = numeric(),
            NES = numeric(), leadingEdge = list())
    }, .package = "fgsea")
    ordered <- list(genes = genes, metabolites = metabolites)
    for (tables in list(ordered, rev(ordered))) {
        findRegionalPathways(object, ident = "region", DE.list = tables,
            database = db, gene_mapping = "ramp", organism = "custom",
            min_path_size = 1, verbose = FALSE)
    }
    half <- length(seen) / 2
    expect_gt(half, 0)
    expect_identical(seen[seq_len(half)], seen[half + seq_len(half)])
    expect_true(any(startsWith(names(seen[[1]]), "RAMP_C_")))
    expect_true(any(startsWith(names(seen[[1]]), "RAMP_G_")))
})
