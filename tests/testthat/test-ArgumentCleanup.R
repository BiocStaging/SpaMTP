test_that("technical fold-change ratios reject invalid values before pooling", {
    data(spamtp_example)
    object <- spamtp_example$experiment
    for (value in list(0.5, 0, -1, NA_real_, Inf, NaN, numeric(), c(1, 2), "2", TRUE, 1 + 1i)) {
        expect_error(findAllDEMs(object, "region", method = "technical",
            foldChangeThreshold = value, verbose = FALSE), "fold-change ratio >= 1")
        expect_error(runDE(NULL, NULL, foldChangeThreshold = value),
            "fold-change ratio >= 1")
    }
    withr::local_seed(901)
    x <- matrix(stats::rpois(40 * 600, 100), 40,
        dimnames = list(paste0("f", seq_len(40)), paste0("p", seq_len(600))))
    x[1:4, 1:300] <- 2 * x[1:4, 1:300]
    object <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = x),
        colData = S4Vectors::DataFrame(region = factor(rep(c("high", "low"), each = 300))))
    for (threshold in c(1, 1.2)) {
        expect_warning(result <- findAllDEMs(object, "region", method = "technical",
            foldChangeThreshold = threshold, verbose = FALSE), "technical subdivisions")
        hits <- result$DEMs[result$DEMs$cluster == "high" &
            result$DEMs$gene %in% paste0("f", 1:4), ]
        expect_equal(nrow(hits), 4L)
        expect_true(all(hits$FDR < 0.05 & hits$logFC > 0.9 & hits$logFC < 1.1))
    }
})

test_that("feature summaries combine paired modalities using their aligned values", {
    data(spamtp_example)
    object <- spamtp_example$experiment
    rna <- SingleCellExperiment::altExp(object, "transcriptome")
    features <- c(rownames(object)[1], rownames(rna)[1:2])
    observed <- summarizeFeatures(object, features, groupBy = "region")
    for (feature in features) {
        source <- if (feature %in% rownames(object)) object else rna
        values <- as.numeric(SummarizedExperiment::assay(source)[feature, ])
        expected <- tapply(values, object$region, mean)
        actual <- observed[observed$feature == feature, ]
        expect_equal(actual$mean, as.numeric(expected[actual$group]))
    }
    expect_error(summarizeFeatures(object, features, groupBy = "region", assay = "main"),
        "existing rows")
    SingleCellExperiment::altExp(object, "duplicate") <- rna
    expect_error(summarizeFeatures(object, features, groupBy = "region"),
        "several alternative experiments")
})

test_that("index ion accessors preserve mass and candidate identities", {
    data(spamtp_example)
    index <- buildMZAnnotationIndex(spamtp_example$metabolites,
        adducts = c("M+H", "M+Na"))
    ions <- annotationIndexIons(index)
    expect_true(all(diff(ions$expected_mz) >= 0))
    expect_setequal(ions$adduct, c("M+H", "M+Na"))
    candidates <- queryMZAnnotationIndex(ions$expected_mz, index, ppm = 0,
        checkFormulaMass = FALSE, checkIsotopes = FALSE, checkAdductNetwork = FALSE)
    expected <- paste(ions$formula, ions$adduct, ions$expected_mz)
    actual <- paste(candidates$formula, candidates$adduct, candidates$expected_mz)
    expect_setequal(actual, expected)
    expect_error(annotationIndexIons(list()), "buildMZAnnotationIndex")
})

test_that("workflow metrics expose the computed representation comparisons", {
    data(spamtp_example)
    fit <- runSpaMTPWorkflow(spamtp_example$experiment,
        modalities = list(main = list(type = "metabolomics", layer = "counts")),
        group = "region", npcs = 2, clusters = 2, associationFeatures = 0,
        structure = list(spatial = FALSE, umap = FALSE, stability_blocks = 0))
    metrics <- workflowMetrics(fit)
    expect_s3_class(metrics, "data.frame")
    expect_true(all(c("representation", "k", "approx_silhouette") %in% names(metrics)))
    expect_true(any(metrics$representation == "main / PCA" & metrics$k == 2))
    expect_error(workflowMetrics(list()), "computed structure metrics")
})
