current_annotation_fixture <- function() {
    data.frame(
        observed_mz = 149.04445,
        Adduct = "M+H",
        Ramp_IDs = "RAMP_C_000219574",
        IsomerNames = "legacy display name",
        Score = 0.95,
        MassScore = 0.99,
        ChemicalScore = 1,
        IsotopeScore = NA_real_,
        AdductNetworkScore = NA_real_,
        Error = 0.2,
        stringsAsFactors = FALSE
    )
}

annotation_object_fixture <- function() {
    counts <- Matrix::Matrix(
        matrix(
            1:4,
            nrow = 2,
            dimnames = list(c("mz-149.04445", "mz-150"), c("spot1", "spot2"))
        ),
        sparse = TRUE
    )
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(mz = c(149.04445, 150)),
        spatialCoords = cbind(x = 1:2, y = c(0, 1))
    )
    SingleCellExperiment::mainExpName(object) <- "RNA"
    object
}

test_that("current annotation store is preferred over legacy db_3", {
    current <- current_annotation_fixture()
    legacy <- data.frame(
        observed_mz = 149.04445,
        Isomers_IDs = "hmdb:HMDB0000606",
        stringsAsFactors = FALSE
    )
    tools <- list(
        mz_annotation = SpaMTP:::.annotation_store(current),
        db_3 = legacy
    )

    result <- SpaMTP:::.select_mz_annotations(tools, "current")

    expect_equal(result$Ramp_IDs, "RAMP_C_000219574")
    expect_equal(attr(result, "annotation_source"), "mz_annotation")
    expect_equal(
        attr(result, "annotation_metadata")$engine,
        "indexed-chemical-v2"
    )
})

test_that("current mode refuses an unversioned legacy annotation", {
    tools <- list(db_3 = data.frame(
        observed_mz = 149.04445,
        Isomers_IDs = "hmdb:HMDB0000606",
        stringsAsFactors = FALSE
    ))

    expect_error(
        SpaMTP:::.select_mz_annotations(tools, "current"),
        "Run annotateSM"
    )
    expect_warning(
        legacy <- SpaMTP:::.select_mz_annotations(tools, "auto"),
        "legacy"
    )
    expect_equal(attr(legacy, "annotation_source"), "db_3-legacy")
})

test_that("current Ramp IDs map directly to current chemical metadata", {
    object <- annotation_object_fixture()
    object <- SpaMTP:::.setStoredData(
        object,
        "mz_annotation",
        SpaMTP:::.annotation_store(current_annotation_fixture())
    )
    chemical <- data.frame(
        ramp_id = "RAMP_C_000219574",
        chem_source_id = "hmdb:HMDB0000606",
        common_name = "D-2-Hydroxyglutaric acid",
        stringsAsFactors = FALSE
    )

    result <- SpaMTP:::.resolve_pathway_metabolite_annotations(
        object,
        annotationSource = "current", chemicalProperties = chemical
    )

    expect_equal(result$ramp_id, "RAMP_C_000219574")
    expect_equal(result$chem_source_id, "hmdb:HMDB0000606")
    expect_equal(result$IsomerNames, "D-2-Hydroxyglutaric acid")
    expect_equal(result$Score, 0.95)
})

test_that("pathway resolver keeps the best scored adduct hypothesis", {
    object <- annotation_object_fixture()
    annotations <- current_annotation_fixture()[c(1, 1), ]
    annotations$Adduct <- c("M+Na", "M+H")
    annotations$Score <- c(0.2, 0.9)
    annotations$Error <- c(4.5, 0.4)
    object <- SpaMTP:::.setStoredData(
        object, "mz_annotation", SpaMTP:::.annotation_store(annotations)
    )
    chemical <- data.frame(
        ramp_id = "RAMP_C_000219574",
        chem_source_id = "hmdb:HMDB0000606",
        common_name = "D-2-Hydroxyglutaric acid",
        stringsAsFactors = FALSE
    )

    result <- SpaMTP:::.resolve_pathway_metabolite_annotations(
        object,
        annotationSource = "current", chemicalProperties = chemical
    )

    expect_equal(nrow(result), 1L)
    expect_equal(result$Adduct, "M+H")
    expect_equal(result$Score, 0.9)
})

test_that("pathway annotation score threshold is user adjustable", {
    object <- annotation_object_fixture()
    annotations <- current_annotation_fixture()[c(1, 1), ]
    annotations$observed_mz <- c(149.04445, 150.00000)
    annotations$Ramp_IDs <- c("RAMP_C_HIGH", "RAMP_C_LOW")
    annotations$Score <- c(0.8, 0.02)
    object <- SpaMTP:::.setStoredData(
        object, "mz_annotation", SpaMTP:::.annotation_store(annotations)
    )
    chemical <- data.frame(
        ramp_id = c("RAMP_C_HIGH", "RAMP_C_LOW"),
        chem_source_id = c("chebi:1", "chebi:2"),
        common_name = c("high score", "low score"),
        stringsAsFactors = FALSE
    )

    strict <- SpaMTP:::.resolve_pathway_metabolite_annotations(
        object,
        scoreThreshold = 0.05, chemicalProperties = chemical
    )
    permissive <- SpaMTP:::.resolve_pathway_metabolite_annotations(
        object,
        scoreThreshold = 0.01, chemicalProperties = chemical
    )

    expect_equal(strict$ramp_id, "RAMP_C_HIGH")
    expect_setequal(permissive$ramp_id, c("RAMP_C_HIGH", "RAMP_C_LOW"))
    expect_equal(
        attr(strict, "annotation_metadata")$analysis_score_threshold,
        0.05
    )
})

test_that("annotationInfo reports legacy objects without treating them as current", {
    object <- annotation_object_fixture()
    object <- SpaMTP:::.setStoredData(
        object,
        "db_3",
        data.frame(
            observed_mz = 100,
            Isomers_IDs = "chebi:1",
            stringsAsFactors = FALSE
        )
    )

    info <- annotationInfo(object)

    expect_equal(info$schema_version, 1L)
    expect_equal(info$engine, "legacy-mass-match")
})

test_that("curated FMP10 annotations are mapped into the current RaMP schema", {
    fmp <- data.frame(
        mass = c(179.0831, 200),
        observed_mz = c(179.0832, 200.001),
        Adduct = c("[M+K]", "[M+FMP10]+"),
        Formula = c("C9H16O", "C5H10O5"),
        Isomers = c("LMFA05000118", "custom"),
        Isomers_IDs = c("LIPIDMAPS:LMFA05000118", "HMDB:HMDBTEST"),
        IsomerNames = c("test lipid", "test metabolite"),
        stringsAsFactors = FALSE
    )
    chemical <- data.frame(
        ramp_id = c("RAMP_C_000000123", "RAMP_C_000000456"),
        chem_source_id = c("lipidmaps:LMFA05000118", "hmdb:HMDBTEST"),
        stringsAsFactors = FALSE
    )

    result <- SpaMTP:::.fmp10_current_annotations(
        fmp,
        massThreshold = 0.05,
        chemicalProperties = chemical
    )
    store <- SpaMTP:::.annotation_store(
        result,
        metadata = list(engine = "curated-fmp10-v2")
    )

    expect_true(SpaMTP:::.annotation_has_current_schema(result))
    expect_setequal(
        result$Ramp_IDs,
        c("RAMP_C_000000123", "RAMP_C_000000456")
    )
    expect_true(all(result$Score > 0 & result$Score <= 1))
    expect_equal(store$metadata$engine, "curated-fmp10-v2")
})

test_that("annotation search falls back to current stored labels", {
    object <- annotation_object_fixture()
    object <- SpaMTP:::.setStoredData(
        object,
        "mz_annotation",
        SpaMTP:::.annotation_store(current_annotation_fixture())
    )

    result <- searchAnnotations(object, "display name", assay = "RNA")

    expect_equal(nrow(result), 1L)
    expect_equal(result$mz_names, "mz-149.04445")
})

test_that("FMP10 exact ties keep names, identities and scores together", {
    object <- annotation_object_fixture()
    SummarizedExperiment::rowData(object)$mz <- c(149.04445, 149.04445)
    reference <- data.frame(
        mass = 149.04445, annotation = "test derivative", Adduct = "M+FMP10",
        Formula = "C5H8O5", Isomers = "test", Isomers_IDs = "HMDB:TEST",
        IsomerNames = "test derivative"
    )
    database <- list(filtered_fmp10 = reference, chem_props = data.frame(
        ramp_id = "RAMP_C_TEST", chem_source_id = "hmdb:TEST"
    ))
    S4Vectors::metadata(object)$modality_annotations <- list(main = list(
        mz_annotation = SpaMTP:::.annotation_store(current_annotation_fixture())
    ))
    annotated <- addFMP10Annotations(object, database = database)
    candidates <- annotationCandidates(annotated)
    expect_identical(candidates$feature, rownames(object))
    expect_identical(length(unique(candidates$feature_id)), 2L)
    expect_identical(length(unique(candidates$ID)), 2L)
    expect_identical(
        summarizeAnnotationScores(annotated, 0)$matched_features, 2L
    )
    metadata <- SummarizedExperiment::rowData(annotated)
    expect_identical(as.character(metadata$all_IsomerNames),
        rep("test derivative", 2)
    )
    expect_identical(as.character(metadata$all_Ramp_IDs),
        rep("RAMP_C_TEST", 2)
    )
    expect_equal(as.numeric(metadata$all_Scores),
        candidates$Score
    )
    expect_identical(annotationCandidates(
        SpaMTP:::.pathwayAnnotationObject(annotated, "main")
    ), candidates)

    renamed <- annotated[2:1, ]
    keys <- SummarizedExperiment::rowData(renamed)$.spamtp_feature_id
    rownames(renamed) <- c("renamed_second", "renamed_first")
    again <- addFMP10Annotations(renamed, database = database)
    expect_identical(
        SummarizedExperiment::rowData(again)$.spamtp_feature_id, keys
    )
    expect_identical(annotationCandidates(again)$feature_id,
        rev(candidates$feature_id)
    )
    expect_equal(annotationCandidates(again)$Score, rev(candidates$Score))
    expect_identical(nrow(annotationCandidates(again[1, ])), 1L)
    expect_identical(nrow(addFMP10Annotations(object,
        database = database, returnOnlyAnnotated = TRUE
    )), 2L)
})

test_that("custom mass matching keeps the nearest rule outside exact ties", {
    object <- annotation_object_fixture()
    custom <- data.frame(annotation = c("near first", "near second"),
        mass = c(149.04446, 149.99999)
    )
    result <- addCustomMZAnnotations(object, custom)
    metadata <- SummarizedExperiment::rowData(result)
    expect_identical(as.character(metadata$all_IsomerNames),
        custom$annotation
    )
})

test_that("FMP10 alternative experiments preserve the primary legacy store", {
    object <- annotation_object_fixture()
    legacy <- data.frame(observed_mz = 149.04445, Isomers_IDs = "HMDB:OLD")
    S4Vectors::metadata(object)$db_3 <- legacy
    SingleCellExperiment::altExp(object, "fmp") <- annotation_object_fixture()
    reference <- data.frame(
        mass = 149.04445, annotation = "test derivative", Adduct = "M+FMP10",
        Formula = "C5H8O5", Isomers = "test", Isomers_IDs = "HMDB:TEST",
        IsomerNames = "test derivative"
    )
    database <- list(filtered_fmp10 = reference, chem_props = data.frame(
        ramp_id = "RAMP_C_TEST", chem_source_id = "hmdb:TEST"
    ))
    annotated <- addFMP10Annotations(object, assay = "fmp", database = database)
    primary <- SpaMTP:::.pathwayAnnotationObject(annotated, "main")
    alternative <- SpaMTP:::.pathwayAnnotationObject(annotated, "fmp")
    expect_identical(S4Vectors::metadata(primary)$db_3, legacy)
    expect_null(S4Vectors::metadata(primary)$mz_annotation)
    expect_identical(annotationCandidates(alternative)$feature, "mz-149.04445")
    expect_identical(annotationCandidates(alternative)$Ramp_IDs, "RAMP_C_TEST")
    expect_identical(summarizeAnnotationScores(alternative, 0)$matched_features,
        1L
    )
})
