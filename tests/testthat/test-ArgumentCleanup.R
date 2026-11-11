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

test_that("functions and S4 methods reject misspelt or legacy arguments in ...", {
    data(spamtp_example)
    spe <- spamtp_example$experiment
    db <- spamtp_example$database

    expect_error(
        fishersPathwayAnalysis(
            list(genes = c("GENE1", "GENE2")),
            universe = list(genes = paste0("GENE", seq_len(8))),
            database = db, minPathSize = 2, pval_cutoff = 0.01, verbose = FALSE
        ),
        "pvalCutoff"
    )
    expect_error(
        fishersPathwayAnalysis(
            list(genes = c("GENE1", "GENE2")),
            universe = list(genes = paste0("GENE", seq_len(8))),
            database = db, minPathSize = 2, ppm = 5, ppm = 10, verbose = FALSE
        ),
        "Duplicate argument"
    )
    expect_error(
        asSpatialExperiment(spe, "unnamed"),
        "Unused unnamed"
    )

    empty_df <- data.frame(analytes_in_pathways = c(1L, 2L), p_val = c(0.01, 0.02))
    expect_error(
        visualisePathways(NULL, empty_df, minN = 3, min_n = 1, verbose = FALSE),
        "minN"
    )
    single_df <- data.frame(
        pathway_name = "Single pathway",
        pathway_id = "WP1",
        p_val = 0.01,
        fdr = 0.02,
        analytes_in_pathways = 5L,
        total_in_pathways = 8L
    )
    expect_error(
        visualisePathways(NULL, single_df, minN = 1, min_n = 1, verbose = FALSE),
        "minN"
    )

    expect_error(
        plotSpatialFeature(spe, "ion1", legend_title = "x"),
        "legendTitle"
    )
    msi <- asCardinal(spe)
    msi_feature <- rownames(asSpatialExperiment(msi))[[1L]]
    expect_error(
        plotSpatialFeature(msi, msi_feature, legend_title = "x"),
        "legendTitle"
    )
    msi_plot <- plotSpatialFeature(
        msi, msi_feature, legendTitle = "Custom title"
    )
    expect_identical(msi_plot$labels$colour, "Custom title")

    expect_error(asSpatialExperiment(spe, assay_name = "counts"), "Unused argument")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(spe),
        colData = cbind(
            SummarizedExperiment::colData(spe),
            SpatialExperiment::spatialCoords(spe)
        )
    )
    expect_error(asSpatialExperiment(se, assay_name = "counts"), "Unused argument")
    expect_error(asCardinal(msi, assay_name = "counts"), "Unused argument")

    rna <- SummarizedExperiment::assay(
        SingleCellExperiment::altExp(spe, "transcriptome")
    )
    expect_error(
        addTranscriptome(spe, rna, name = "extra", assay_name = "counts"),
        "assayName"
    )

    img_file <- tempfile(fileext = ".png")
    on.exit(unlink(img_file), add = TRUE)
    png::writePNG(matrix(0.5, 4, 4), img_file)
    expect_error(
        addSpatialImage(spe, img_file, sample_id = "s1"),
        "sampleId"
    )
    expect_error(
        addSMImage(img_file, spe, interactive = FALSE, scale_factor = 2),
        "scaleFactor"
    )

    precomputed <- multiOmicIntegration(
        spe, dimsList = list(seq_len(2), seq_len(2)), returnIntermediate = TRUE
    )
    expect_error(
        multiOmicIntegration(precomputed, dims.list = list(seq_len(2), seq_len(2))),
        "dimsList"
    )

    idx <- buildPathwayIndex(db)
    expect_error(
        plotPathwaysSpatially(
            "P1", spe, assay = "transcriptome", pathwayIndex = idx,
            min.cutoff = -1
        ),
        "minCutoff"
    )
})
