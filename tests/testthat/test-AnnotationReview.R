# Regression cases from the annotation accuracy and usability review.
review_spe <- function(mz) {
    counts <- matrix(seq_len(length(mz) * 2L), nrow = length(mz),
        dimnames = list(paste0("ion", seq_along(mz)), c("p1", "p2")))
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts), rowData = S4Vectors::DataFrame(mz = mz),
        spatialCoords = cbind(x = 1:2, y = c(0, 0)))
}

review_db <- function(formula = "C6H12O6") {
    data.frame(formula = formula, exactmass = MetaboCoreUtils::calculateMass(formula),
        isomers_names = paste0("candidate ", seq_along(formula)))
}

test_that("accurate CHO acids and large lipids retain their chemical evidence", {
    for (mode in c("positive", "negative")) {
        adduct <- if (mode == "positive") "M+H" else "M-H"
        index <- buildMZAnnotationIndex(review_db(c("C2H2O4", "C60H116O6", "C61H114O6")),
            polarity = mode, adducts = adduct, infer_structure = "never")
        result <- queryMZAnnotationIndex(index$expected_mz, index, ppm = 0)
        expect_equal(result$chemical_score, rep(1, 3))
        expect_equal(result$mass_score, rep(1, 3))
    }
})

test_that("FMP10 net product masses include the observed cation charge", {
    dopamine <- MetaboCoreUtils::calculateMass("C8H11NO2")
    rules <- maldiMatrixRules("FMP-10")
    names <- c("M+FMP10", "M+2FMP10a", "M+2FMP10b")
    actual <- dopamine + rules$mass_shift[match(names, rules$name)]
    # Independently specified product formulas, minus one electron each.
    expected <- MetaboCoreUtils::calculateMass(c("C28H25N2O2", "C47H36N3O2", "C48H38N3O2")) - 0.000548579909065
    expect_equal(actual, unname(expected), tolerance = 1e-8)
    expect_equal(round(actual, 4), c(421.1911, 674.2802, 688.2959))
})

test_that("equal-mass features receive one copy of each candidate", {
    index <- buildMZAnnotationIndex(review_db(), adducts = "M+H")
    mz <- index$expected_mz
    table <- annotateTable(data.frame(row_id = c(8, 2, 5), mz = rep(mz, 3)),
        index = index, ppm_error = 0, verbose = FALSE)
    expect_equal(nrow(table), 3L)
    expect_setequal(table$ID, c(8, 2, 5))
    obj <- annotateSM(review_spe(rep(mz, 3)), index = index, ppm_error = 0,
        return.only.annotated = FALSE, verbose = FALSE)
    expect_equal(nrow(annotationCandidates(obj)), 3L)
    expect_equal(annotationInfo(obj)$candidates, 3L)
    expect_identical(as.character(SummarizedExperiment::rowData(obj)$all_IsomerNames), rep("candidate 1", 3))
    big <- annotateBigData(rep(mz, 3), index = index, ppm_error = 0, verbose = FALSE)
    expect_equal(nrow(big), 1L)
    expect_identical(big$all_IsomerNames, "candidate 1")
})

test_that("index reuse honors restrictions and rejects conflicting construction options", {
    db <- review_db()
    index <- buildMZAnnotationIndex(db, adducts = c("M+H", "M+Na"), infer_structure = "never")
    mz <- index$expected_mz
    hit <- annotateMZ(mz, index = index, adducts = "M+Na", ppm = 0)
    expect_identical(hit$adduct, "M+Na")
    hit <- annotateTable(data.frame(row_id = 1:2, mz = mz), index = index,
        adducts = "[M+Na]+", ppm_error = 0, verbose = FALSE)
    expect_identical(hit$Adduct, "M+Na")
    expect_error(annotateMZ(mz, index = index, adducts = "M+K"), "unavailable")
    expect_error(annotateMZ(mz, index = index, db = db), "omit db")
    expect_error(annotateMZ(mz, index = index, rules = adductRules()), "omit db")
    expect_error(annotateMZ(mz, index = index, infer_structure = "always"), "infer_structure")
    expect_error(annotateMZ(mz, index = index, min_structure_score = 0.9), "min_structure_score")
    expect_silent(annotateMZ(mz, index = index, infer_structure = "never"))
    expect_length(index$expected_mz, 2L)
})

test_that("annotation storage and inspection never fetch database provenance", {
    local_mocked_bindings(.spamtp_db_resource = function(...) stop("unexpected resource fetch"))
    db <- review_db()
    index <- buildMZAnnotationIndex(db, adducts = "M+H", infer_structure = "never")
    obj <- annotateSM(review_spe(index$expected_mz), index = index, ppm_error = 0,
        verbose = FALSE, return.only.annotated = FALSE)
    expect_true(is.na(annotationInfo(obj)$ramp_version))
    S4Vectors::metadata(obj)$mz_annotation <- NULL
    expect_true(is.na(annotationInfo(obj)$ramp_version))
    attr(db, "spamtp_database") <- list(version = "test-release", source = "spamtpdb")
    index <- buildMZAnnotationIndex(db, adducts = "M+H", infer_structure = "never")
    obj <- annotateSM(review_spe(index$expected_mz), index = index, verbose = FALSE)
    expect_identical(annotationInfo(obj)$ramp_version, "test-release")
})

test_that("custom labels append without erasing unmatched annotations", {
    obj <- review_spe(c(100, 200, 500))
    SummarizedExperiment::rowData(obj)$all_IsomerNames <- c("2-HG", "hexose; inositol", "No Annotation")
    SummarizedExperiment::rowData(obj)$all_Scores <- c("0.9", "0.8; 0.8", "No Annotation")
    labels <- data.frame(annotation = c("manual", "candidate A"), mass = c(100, 500.125))
    out <- addCustomMZAnnotations(obj, labels, mass.threshold = 0.125)
    rd <- SummarizedExperiment::rowData(out)
    expect_identical(as.character(rd$all_IsomerNames), c("2-HG; manual", "hexose; inositol", "candidate A"))
    expect_identical(rd$all_Scores, SummarizedExperiment::rowData(obj)$all_Scores)
    expect_equal(as.numeric(rd$mass_diff_da[3]), 0.125)
    expect_equal(as.numeric(rd$ppm_diff[3]), 0.125 / 500.125 * 1e6)
    out <- addCustomMZAnnotations(obj, labels, mass.threshold = 0, append = FALSE)
    expect_identical(as.character(SummarizedExperiment::rowData(out)$all_IsomerNames),
        c("manual", "hexose; inositol", "No Annotation"))
    again <- addCustomMZAnnotations(out, labels[1, ], mass.threshold = 0)
    expect_identical(SummarizedExperiment::rowData(again)$all_IsomerNames,
        SummarizedExperiment::rowData(out)$all_IsomerNames)
})

test_that("curated FMP10 stores use exact reference row identities and ppm errors", {
    local_mocked_bindings(.spamtp_db_resource = function(...) stop("unexpected resource fetch"))
    mass <- 153.0789786 + 268.1126
    reference <- data.frame(mass = mass, annotation = "dopamine derivative", Adduct = "M+FMP10",
        Formula = "C8H11NO2", Isomers = "dopamine", Isomers_IDs = "hmdb:1", IsomerNames = "dopamine derivative")
    database <- list(filtered_fmp10 = reference, chem_props = data.frame(ramp_id = "RAMP_C_1", chem_source_id = "hmdb:1"))
    obj <- review_spe(c(mass, mass, 900))
    out <- addFMP10Annotations(obj, database = database, mass.threshold = 0)
    candidates <- annotationCandidates(out)
    expect_equal(nrow(candidates), 2L)
    expect_identical(candidates$Reference_mz, rep(mass, 2))
    expect_equal(candidates$Score, c(1, 1))
    expect_identical(as.character(SummarizedExperiment::rowData(out)$all_IsomerNames),
        c("dopamine derivative", "dopamine derivative", "No Annotation"))
    obj <- review_spe(mass + 0.001)
    out <- addFMP10Annotations(obj, database = database)
    expect_equal(as.numeric(SummarizedExperiment::rowData(out)$all_Errors),
        round(annotationCandidates(out)$Error, 4))
    empty <- addFMP10Annotations(review_spe(900), database = database)
    expect_equal(nrow(annotationCandidates(empty)), 0L)
})

test_that("single-pixel profiles report unavailable correlations without crashing", {
    isotopes <- flagIsotopePeaks(c(399.1436, 400.1466, 410.2), matrix(c(100, 14, 50), ncol = 1))
    groups <- groupAdductPeaks(c(149.0444, 171.0264, 300.1), matrix(c(5, 3, 1), ncol = 1))
    expect_equal(nrow(isotopes), 3L)
    expect_false(any(isotopes$isotope))
    expect_true(is.list(groups))
})

test_that("literal annotation lookup and metadata helpers handle edge cases", {
    obj <- review_spe(c(100, 200, 300, 400, 500))
    SummarizedExperiment::rowData(obj)$all_IsomerNames <- c("2-{[amino]}acetic acid; a|b", "banana; 17xbeta", "17.beta", "No Annotation", "No Annotation")
    expect_identical(rownames(searchAnnotations(obj, "a|b")), "ion1")
    expect_identical(rownames(searchAnnotations(obj, "2-{[amino]}acetic acid", search.exact = TRUE)), "ion1")
    expect_identical(rownames(searchAnnotations(obj, "17.beta")), "ion3")
    expect_identical(findDuplicateAnnotations(obj), character())
    expect_error(getMZMetadata(obj, "missing"), "Unknown m/z feature: 'missing'")
    expect_error(getMZMetadata(obj, "ion1", metadata.column = "missing"), "metadata column")
    expect_length(getMZMetadata(obj, "ion1"), 2L)
    expect_equal(nrow(getRefinedAnnotations(obj, n = 0)), 0L)
    expect_identical(rownames(getRefinedAnnotations(obj, n = 1)), "ion3")
    expect_error(getRefinedAnnotations(obj, n = -1), "non-negative integer")
})

test_that("refinement counts stored candidates rather than semicolons in names", {
    db <- review_db()
    db$isomers_names <- "hexose; glucose"
    index <- buildMZAnnotationIndex(db, adducts = "M+H")
    obj <- annotateSM(review_spe(c(index$expected_mz, 999)), index = index,
        return.only.annotated = FALSE, verbose = FALSE)
    expect_identical(rownames(getRefinedAnnotations(obj, n = 1)), "ion1")
})

test_that("functional group perception separates amines from nitro and sulfonamide groups", {
    result <- deconvolveSMILES(c("c1ccccc1[N+](=O)[O-]", "NS(=O)(=O)c1ccccc1", "C[N+](C)(C)CCO", "OP(=O)(O)O", "NCCc1ccc(O)c(O)c1"))
    expect_equal(result$primary_amine_sites, c(0, 0, 0, 0, 1))
    expect_equal(result$tertiary_amine_sites, rep(0, 5))
    expect_equal(result$basic_sites[1:4], rep(0, 4))
    expect_equal(result$fmp10_reactive_sites, c(0, 0, 0, 0, 3))
    expect_equal(result$alcohol_hydroxyl_sites[4], 0)
    expect_equal(result$phosphate_acid_sites[4], 3)
    expect_equal(result$positive_mode_score[3], 1)
})

test_that("isotope evidence accounts for multimers and derivatization carbons", {
    for (adduct in c("2M+H", "M+ACN+H")) {
        index <- buildMZAnnotationIndex(review_db(), adducts = adduct)
        carbon <- if (adduct == "2M+H") 12 else 8
        mz <- index$expected_mz
        spectrum <- data.frame(mz = c(mz, mz + 1.00335483507), intensity = c(100, 100 * carbon * 0.0107))
        out <- queryMZAnnotationIndex(mz, index, ms1_spectrum = spectrum, check_adduct_network = FALSE)
        expect_equal(out$isotope_score, 1)
    }
    db <- review_db("C8H11NO2")
    db$fmp10_reactive_sites <- 3
    index <- buildMZAnnotationIndex(db, maldi_matrix = "FMP-10", adducts = "M+FMP10")
    mz <- index$expected_mz
    out <- queryMZAnnotationIndex(mz, index,
        ms1_spectrum = data.frame(mz = c(mz, mz + 1.00335483507), intensity = c(100, 29.96)))
    expect_equal(out$isotope_score, 1)
})

test_that("adduct aliases are deduplicated and FDR supports numeric spectra", {
    index <- buildMZAnnotationIndex(review_db(), adducts = c("M+H", "[M+H]+"))
    expect_length(index$expected_mz, 1L)
    mz <- index$expected_mz
    expect_equal(nrow(queryMZAnnotationIndex(mz, index)), 1L)
    expect_error(estimateAnnotationFDR(mz, index, ppm = 0), "Supply non-zero shifts_da")
    expect_s3_class(estimateAnnotationFDR(mz, index, ppm = 0, shifts_da = c(-1, 1)), "spamtp_annotation_fdr")
    a <- estimateAnnotationFDR(mz, index, ms1_spectrum = c(mz + 1, mz))
    b <- estimateAnnotationFDR(mz, index, ms1_spectrum = data.frame(mz = c(mz + 1, mz), intensity = 1))
    expect_equal(a, b)
})

test_that("annotation rejects mistyped mass columns and missing feature names", {
    index <- buildMZAnnotationIndex(review_db(), adducts = "M+H")
    obj <- review_spe(index$expected_mz)
    expect_error(annotateSM(obj, index = index, raw.mz.column = "typo"), "column 'typo'")
    rownames(obj) <- NULL
    expect_error(annotateSM(obj, index = index), "feature rownames")
})

test_that("explicit structure parsing choices do not fetch precomputed tables", {
    calls <- 0L
    local_mocked_bindings(.spamtp_db_resource = function(...) {
        calls <<- calls + 1L
        stop("unexpected resource fetch")
    })
    db <- review_db("C6H7NO2S")
    db$smiles <- "NS(=O)(=O)c1ccccc1"
    attr(db, "spamtp_database") <- list(source = "spamtpdb", version = "old")
    buildMZAnnotationIndex(db, adducts = "M+H", infer_structure = "never")
    db$structure_valid <- TRUE
    db$primary_amine_sites <- 1
    db$fmp10_reactive_sites <- 1
    refreshed <- buildMZAnnotationIndex(db, adducts = "M+H", infer_structure = "always")
    expect_equal(calls, 0L)
    expect_equal(refreshed$compounds$primary_amine_sites, 0)
    expect_equal(refreshed$compounds$fmp10_reactive_sites, 0)
})
