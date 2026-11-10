# Development assistance: see inst/COPYRIGHTS.

annotatedFixture <- function() {
    chemical <- data.frame(
        ramp_id = paste0("RAMP_C_", 1:3), chem_source_id = paste0("hmdb:HMDB", 1:3),
        common_name = c("compoundA", "compoundB", "compoundC"), chem_data_source = "hmdb",
        mol_formula = c("C5H8O5", "C5H8O5", "C6H12O6"),
        monoisotop_mass = c(148.037173366, 148.037173366, 180.063388104)
    )
    index <- buildMZAnnotationIndex(chemical, adducts = "M+H", infer_structure = "never")
    object <- nativeFixture()
    SummarizedExperiment::rowData(object)$mz <- c(
        MetaboCoreUtils::mass2mz(c(148.037173366, 180.063388104), "[M+H]+")[, 1], 600
    )
    annotateSM(object,
        index = index, ppm_error = 5, min_score = 0,
        return.only.annotated = FALSE, save.intermediate = TRUE, verbose = FALSE
    )
}

test_that("findNearestMZ matches each target mass", {
    object <- nativeFixture()
    expect_identical(findNearestMZ(object, c(201, 99, 290)), c("b", "a", "constant"))
    expect_error(findNearestMZ(object, NA_real_), "finite")
    expect_error(findNearestMZ(object, "100"), "finite")
})

test_that("annotation candidates are ranked per feature with identity counts", {
    annotated <- annotatedFixture()
    candidates <- annotationCandidates(annotated)
    expect_identical(unique(candidates$feature), c("a", "b"))
    expect_true(all(candidates$rank == stats::ave(seq_len(nrow(candidates)),
        candidates$feature,
        FUN = seq_along
    )))
    expect_identical(unique(candidates$identities[candidates$feature == "a"]), 2L)
    expect_identical(unique(candidates$identities[candidates$feature == "b"]), 1L)
    expect_true(all(diff(candidates$Score[candidates$feature == "a"]) <= 0))
    expect_identical(unique(annotationCandidates(annotated, features = 181.07)$feature), "b")
    expect_identical(nrow(annotationCandidates(annotated, features = "constant")), 0L)
    expect_identical(nrow(annotationCandidates(annotated, minScore = 2)), 0L)
    expect_error(annotationCandidates(annotated, features = "absent"), "Unknown features")

    summary <- summarizeAnnotationScores(annotated, thresholds = c(0, 2))
    expect_identical(summary$matched_features, c(2L, 0L))
    expect_identical(summary$single_identity_features, c(1L, 0L))
    expect_identical(summary$ambiguous_features, c(1L, 0L))
    expect_s3_class(plotAnnotationScores(annotated, threshold = 0.05), "ggplot")
})

test_that("candidate features follow the stored masses after subsetting", {
    annotated <- annotatedFixture()
    reordered <- annotated[c("b", "constant", "a"), ]
    candidates <- annotationCandidates(reordered)
    expect_identical(unique(candidates$feature), c("b", "a"))
    expect_equal(
        candidates$observed_mz,
        SummarizedExperiment::rowData(reordered)[candidates$feature, "mz"]
    )
})

test_that("candidate tables, summaries and plots exclude removed features", {
    annotated <- annotatedFixture()
    original <- annotationCandidates(annotated)
    retained <- annotated["b", , drop = FALSE]
    candidates <- annotationCandidates(retained)
    expect_identical(unique(candidates$feature), "b")
    expect_equal(candidates$Score, original$Score[original$feature == "b"])

    summary <- summarizeAnnotationScores(retained, thresholds = c(0, 2))
    expect_identical(summary$candidates, c(nrow(candidates), 0L))
    expect_identical(summary$matched_features, c(1L, 0L))
    expect_identical(summary$single_identity_features, c(1L, 0L))
    expect_identical(summary$ambiguous_features, c(0L, 0L))
    expect_equal(plotAnnotationScores(retained)$data$score, max(candidates$Score))
    expect_identical(annotationCandidates(annotated), original)

    unmatched <- annotated["constant", , drop = FALSE]
    empty <- annotated[FALSE, , drop = FALSE]
    expect_equal(nrow(annotationCandidates(unmatched)), 0L)
    expect_equal(nrow(annotationCandidates(empty)), 0L)
    expect_equal(summarizeAnnotationScores(unmatched)$matched_features, rep(0L, 4))
    expect_equal(summarizeAnnotationScores(empty)$matched_features, rep(0L, 4))
    expect_equal(nrow(plotAnnotationScores(unmatched)$data), 0L)
    expect_equal(nrow(plotAnnotationScores(empty)$data), 0L)
})

test_that("older candidate stores can be filtered by the retained masses", {
    annotated <- annotatedFixture()
    metadata <- S4Vectors::metadata(annotated)
    metadata$mz_annotation$results$feature <- NULL
    metadata$mz_annotation$results$feature_id <- NULL
    metadata$db_3$feature <- NULL
    metadata$db_3$feature_id <- NULL
    S4Vectors::metadata(annotated) <- metadata
    retained <- annotated[c("constant", "b"), , drop = FALSE]
    expect_identical(unique(annotationCandidates(retained)$feature), "b")
    expect_identical(summarizeAnnotationScores(retained, 0)$matched_features, 1L)
    expect_equal(nrow(plotAnnotationScores(retained)$data), 1L)
    expect_equal(nrow(annotationCandidates(annotated["constant", ])), 0L)
})

duplicateMassFixture <- function() {
    object <- annotatedFixture()
    SummarizedExperiment::rowData(object)$mz[2] <-
        SummarizedExperiment::rowData(object)$mz[1]
    object <- annotateSM(object,
        db = data.frame(
            ramp_id = "RAMP_C_1", chem_source_id = "hmdb:HMDB1",
            common_name = "compoundA", chem_data_source = "hmdb",
            mol_formula = "C5H8O5", monoisotop_mass = 148.037173366
        ),
        raw.mz.column = "mz", adducts = "M+H", ppm_error = 5,
        min_score = 0, return.only.annotated = FALSE, verbose = FALSE,
        infer_structure = "never"
    )
    object
}

test_that("stored feature keys distinguish duplicate measured masses", {
    object <- duplicateMassFixture()
    original <- annotationCandidates(object)
    retained <- annotationCandidates(object["b", , drop = FALSE])
    expect_identical(unique(retained$feature), "b")
    expect_equal(retained$ID, original$ID[original$feature == "b"])
    expect_equal(nrow(retained), sum(original$feature == "b"))
})

test_that("annotation identity survives renaming and custom annotation", {
    object <- annotatedFixture()
    original <- annotationCandidates(object)
    keys <- SummarizedExperiment::rowData(object)$.spamtp_feature_id
    rownames(object) <- c("b", "a", "renamed_constant")
    candidates <- annotationCandidates(object)
    expect_identical(candidates$feature,
        unname(c(a = "b", b = "a")[original$feature])
    )
    expect_equal(candidates$observed_mz, original$observed_mz)
    expect_identical(summarizeAnnotationScores(object, 0)$matched_features, 2L)
    expect_equal(sort(plotAnnotationScores(object)$data$score),
        sort(plotAnnotationScores(annotatedFixture())$data$score)
    )

    custom <- data.frame(annotation = "custom", mass = original$observed_mz[1],
        .spamtp_feature_id = "must not overwrite the key"
    )
    object <- addCustomMZAnnotations(object, custom,
        annotation.column = "custom_name"
    )
    expect_identical(
        SummarizedExperiment::rowData(object)$.spamtp_feature_id, keys
    )
    expect_identical(annotationCandidates(object)$feature, candidates$feature)
    retained <- object[c("renamed_constant", "a"), ]
    expect_identical(unique(annotationCandidates(retained)$feature), "a")
    expect_equal(annotationCandidates(retained)$observed_mz,
        original$observed_mz[original$feature == "b"]
    )
    expect_error(addCustomMZAnnotations(object, custom,
        annotation.column = ".spamtp_feature_id"
    ), "must not replace")
})

test_that("missing or duplicated preserved identities are rejected", {
    object <- annotatedFixture()
    SummarizedExperiment::rowData(object)$.spamtp_feature_id <- NULL
    expect_error(annotationCandidates(object),
        "keys must be present and unique"
    )
    object <- annotatedFixture()
    SummarizedExperiment::rowData(object)$.spamtp_feature_id[2] <- "1"
    expect_error(annotationCandidates(object),
        "keys must be present and unique"
    )
    expect_error(addCustomMZAnnotations(object,
        data.frame(annotation = "custom", mass = 100)
    ), "keys must be present and unique")
})

test_that("legacy equal-mass candidates require unambiguous feature identity", {
    object <- duplicateMassFixture()
    metadata <- S4Vectors::metadata(object)
    for (field in c("feature", "feature_id")) {
        metadata$mz_annotation$results[[field]] <- NULL
        metadata$db_3[[field]] <- NULL
    }
    S4Vectors::metadata(object) <- metadata
    # Masses also match the old positions after swapping the identical masses.
    expect_error(annotationCandidates(object[c("b", "a", "constant"), ]),
        "ambiguous feature masses"
    )
    expect_error(summarizeAnnotationScores(object["b", ]),
        "ambiguous feature masses"
    )
    expect_error(plotAnnotationScores(object), "ambiguous feature masses")
    expect_identical(nrow(annotationCandidates(object["constant", ])), 0L)
    expect_identical(nrow(annotationCandidates(object[FALSE, ])), 0L)
})

test_that("unrelated duplicated masses do not block legacy annotation", {
    object <- annotatedFixture()
    metadata <- S4Vectors::metadata(object)
    metadata$mz_annotation$results$feature_id <- NULL
    metadata$db_3$feature_id <- NULL
    S4Vectors::metadata(object) <- metadata
    object <- object[c("a", "constant", "constant"), ]
    rownames(object) <- c("renamed_a", "unmatched1", "unmatched2")
    expect_identical(unique(annotationCandidates(object)$feature), "renamed_a")
    expect_identical(summarizeAnnotationScores(object, 0)$matched_features, 1L)
})

test_that("conflicting legacy mass columns cannot change row identity", {
    object <- annotatedFixture()
    metadata <- S4Vectors::metadata(object)
    metadata$mz_annotation$results$feature_id <- NULL
    metadata$db_3$feature_id <- NULL
    S4Vectors::metadata(object) <- metadata
    SummarizedExperiment::rowData(object)$raw_mz <-
        SummarizedExperiment::rowData(object)$mz[c(2, 1, 3)]
    expect_error(annotationCandidates(object), "ambiguous feature masses")
})

test_that("duplicate display names cannot merge different annotated rows", {
    object <- annotatedFixture()
    rownames(object)[1:2] <- "same_label"
    expect_error(annotationCandidates(object),
        "unique, non-empty feature names"
    )
    expect_error(summarizeAnnotationScores(object),
        "unique, non-empty feature names"
    )
    expect_error(plotAnnotationScores(object),
        "unique, non-empty feature names"
    )
})

test_that("candidate helpers cannot borrow an alternative experiment's store", {
    object <- annotatedFixture()
    expected <- annotationCandidates(object)
    alternative <- SingleCellExperiment::altExp(object, "transcriptome")
    SummarizedExperiment::rowData(alternative) <-
        SummarizedExperiment::rowData(object)[2:1, ]
    SingleCellExperiment::altExp(object, "transcriptome") <- alternative
    db <- data.frame(
        ramp_id = "RAMP_C_OTHER", chem_source_id = "hmdb:OTHER",
        common_name = "other", chem_data_source = "hmdb",
        mol_formula = "C5H8O5", monoisotop_mass = 148.037173366
    )
    result <- annotateSM(object, assay = "transcriptome", db = db,
        adducts = "M+H", min_score = 0, return.only.annotated = FALSE,
        infer_structure = "never", verbose = FALSE
    )
    expect_identical(annotationCandidates(result), expected)
    expect_identical(summarizeAnnotationScores(result, 0)$matched_features, 2L)
    expect_identical(plotAnnotationScores(result)$data,
        plotAnnotationScores(object)$data
    )
    object <- nativeFixture()
    alternative <- SingleCellExperiment::altExp(object, "transcriptome")
    SummarizedExperiment::rowData(alternative)$mz <- c(149.04445, 200)
    SingleCellExperiment::altExp(object, "transcriptome") <- alternative
    result <- annotateSM(object, assay = "transcriptome", db = db,
        adducts = "M+H", return.only.annotated = FALSE,
        infer_structure = "never", verbose = FALSE
    )
    expect_error(annotationCandidates(result), "No current indexed")
})

test_that("custom labels cannot recreate missing keys against an old store", {
    object <- annotatedFixture()
    SummarizedExperiment::rowData(object)$.spamtp_feature_id <- NULL
    object <- object[c("b", "a", "constant"), ]
    result <- addCustomMZAnnotations(object,
        data.frame(annotation = "custom", mass = 149.04445)
    )
    expect_null(SummarizedExperiment::rowData(result)$.spamtp_feature_id)
    expect_error(annotationCandidates(result), "keys must be present and unique")

    # A rowData-only reannotation must not recreate keys for the old store.
    result <- annotateSM(object, db = data.frame(
        ramp_id = "RAMP_C_OTHER", chem_source_id = "hmdb:OTHER",
        common_name = "other", chem_data_source = "hmdb",
        mol_formula = "C5H8O5", monoisotop_mass = 148.037173366
    ), adducts = "M+H", save.intermediate = FALSE,
    return.only.annotated = FALSE, infer_structure = "never", verbose = FALSE)
    expect_null(SummarizedExperiment::rowData(result)$.spamtp_feature_id)
    expect_error(annotationCandidates(result), "keys must be present and unique")
})

test_that("normalization accepts recorded totals and keeps dense results dense", {
    object <- nativeFixture()
    totals <- c(10, 20, 30, 40, 50, 60)
    object$tic <- totals
    result <- normalizeSMData(object, "LogNormalize", libSizes = "tic", verbose = FALSE)
    counts <- SummarizedExperiment::assay(object, "counts")
    expected <- log1p(sweep(counts, 2, totals, "/") * 1e4)
    expect_equal(as.matrix(SummarizedExperiment::assay(result, "logcounts")), expected,
        ignore_attr = TRUE
    )
    expect_identical(
        SummarizedExperiment::assay(normalizeSMData(object, "RC",
            libSizes = totals,
            verbose = FALSE
        ), "normcounts"),
        SummarizedExperiment::assay(normalizeSMData(object, "RC",
            libSizes = "tic",
            verbose = FALSE
        ), "normcounts")
    )
    expect_error(normalizeSMData(object, libSizes = 1:2), "one finite")
    expect_error(normalizeSMData(object, libSizes = "absent"), "colData")

    dense <- object
    SummarizedExperiment::assay(dense, "counts") <- Matrix::Matrix(counts, sparse = TRUE)
    dense <- normalizeSMData(dense, "LogNormalize", verbose = FALSE)
    expect_true(is.matrix(SummarizedExperiment::assay(dense, "logcounts")))

    sparse <- object
    zeros <- counts
    zeros[, 2:6] <- 0
    SummarizedExperiment::assay(sparse, "counts") <- Matrix::Matrix(zeros, sparse = TRUE)
    sparse <- normalizeSMData(sparse, "LogNormalize", verbose = FALSE)
    expect_s4_class(SummarizedExperiment::assay(sparse, "logcounts"), "sparseMatrix")

    delayed <- object
    SummarizedExperiment::assay(delayed, "counts") <- DelayedArray::DelayedArray(counts)
    delayed <- normalizeSMData(delayed, "LogNormalize", verbose = FALSE)
    expect_s4_class(SummarizedExperiment::assay(delayed, "logcounts"), "DelayedMatrix")
    expect_equal(as.matrix(SummarizedExperiment::assay(delayed, "logcounts")),
        as.matrix(SummarizedExperiment::assay(dense, "logcounts")),
        ignore_attr = TRUE
    )
})

test_that("delayed scaling matches base scale without materialising", {
    object <- nativeFixture()
    values <- SummarizedExperiment::assay(object, "logcounts")
    result <- scaleSMData(object)
    scaled <- SummarizedExperiment::assay(result, "scaled")
    expect_s4_class(scaled, "DelayedMatrix")
    expected <- t(scale(t(values)))
    expected["constant", ] <- 0
    expect_equal(as.matrix(scaled), expected, ignore_attr = TRUE)
    in_memory <- scaleSMData(object, delayed = FALSE)
    expect_true(is.matrix(SummarizedExperiment::assay(in_memory, "scaled")))
    uncentred <- scaleSMData(object, center = FALSE)
    expect_equal(as.matrix(SummarizedExperiment::assay(uncentred, "scaled")),
        t(scale(t(values), center = FALSE)),
        ignore_attr = TRUE
    )
})

test_that("spatial feature plots read paired alternative experiments", {
    object <- nativeFixture()
    plots <- plotSpatialFeature(object, c("a", "gene2"), combine = FALSE)
    expect_length(plots, 2L)
    transcript <- SummarizedExperiment::assay(
        SingleCellExperiment::altExp(object, "transcriptome"), "counts"
    )["gene2", ]
    expect_equal(plots[[2]]$data$intensity, as.numeric(transcript))
    # The primary experiment takes precedence for names present in both.
    expect_equal(
        plots[[1]]$data$intensity,
        as.numeric(SummarizedExperiment::assay(object, "counts")["a", ])
    )
    expect_identical(
        plots[[1]]$labels$colour %||% ggplot2::get_labs(plots[[1]])$colour,
        "counts"
    )
    expect_error(plotSpatialFeature(object, "absent"), "not found")
    expect_error(plotSpatialFeature(object, "a", assayName = "absent"), "not found")
})

test_that("spatial label plots colour and shape pixels by colData", {
    object <- nativeFixture()
    object$side <- rep(c("left", "right"), each = 3)
    plot <- plotSpatialLabels(object, "region",
        shapeBy = "side",
        colours = c(core = "red", edge = "blue")
    )
    expect_s3_class(plot, "ggplot")
    expect_identical(levels(plot$data$label), c("core", "edge"))
    expect_identical(as.character(plot$data$shape), object$side)
    expect_error(plotSpatialLabels(object, "absent"), "colData")
    expect_error(
        plotSpatialLabels(object, "region", colours = c(core = "red")),
        "every label"
    )
})

test_that("marker statistics become complete feature-by-region rank matrices", {
    markers <- list(
        expression = matrix(0, 3, 2, dimnames = list(c("f1", "f2", "f3"), c("A", "B"))),
        DEMs = data.frame(
            gene = c("f1", "f2", "f1", "f3"), cluster = c("A", "A", "B", "B"),
            logFC = c(1.5, -0.5, 0.2, 2), mean_auc = c(0.9, 0.2, 0.6, 0.8)
        )
    )
    ranks <- markerRanks(markers)
    expect_equal(ranks, matrix(c(1.5, -0.5, NA, 0.2, NA, 2), 3,
        dimnames = list(c("f1", "f2", "f3"), c("A", "B"))
    ))
    expect_equal(markerRanks(markers, "mean_auc")["f3", "B"], 0.8)
    expect_error(markerRanks(markers, "absent"), "numeric column")
    expect_error(markerRanks(list()), "findAllDEMs")
})
