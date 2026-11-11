# Development assistance: see inst/COPYRIGHTS.

test_that("complete formulas retain atom counts and incomplete counts stay unknown", {
    expect_equal(.element_count("CH3OH", "H"), 4L)
    expect_equal(.element_count("NaCl", "C"), 0L)
    expect_setequal(.formula_elements("CH3OH"), c("C", "H", "O"))

    formulas <- c("C5H9NO4R2", "C18H28N2O4.HCl", "C6H12O6[59Co]",
        "C6H11[1H]O6", "C20H30D4O5", NA_character_, ""
    )
    for (formula in formulas) {
        expect_silent(counts <- .formula_counts(formula))
        expect_null(counts, info = formula)
        expect_identical(.element_count(formula, "C"), NA_integer_)
    }
})

test_that("formula warnings are summarized once and valid scoring is unchanged", {
    target <- data.frame(
        formula = c("C5H9NO4R2", "C18H28N2O4.HCl", "C6H12O6[59Co]", "C5H8O5"),
        exactmass = c(100, 200, 300, 148.037173366)
    )
    index <- buildMZAnnotationIndex(target,
        adducts = c("M+H", "M+Na"), inferStructure = "never"
    )
    observed <- rep(index$expected_mz, 2)
    spectrum <- data.frame(mz = index$expected_mz, intensity = 100)
    warnings <- character()
    result <- withCallingHandlers(
        queryMZAnnotationIndex(observed, index, ppm = 0,
            ms1Spectrum = spectrum
        ),
        warning = function(condition) {
            warnings <<- c(warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
        }
    )
    expect_length(warnings, 1L)
    expect_match(warnings, "3 distinct matched reference formula")
    expect_equal(nrow(result), length(observed))
    unavailable <- result$formula != "C5H8O5"
    expect_true(all(result$formula_status[unavailable] == "unavailable"))
    expect_true(all(is.na(result$chemical_score[unavailable])))
    expect_true(all(is.na(result$isotope_score[unavailable])))
    expect_true(all(is.finite(result$score)))

    valid_index <- buildMZAnnotationIndex(target[4, ],
        adducts = c("M+H", "M+Na"), inferStructure = "never"
    )
    expect_silent(valid <- queryMZAnnotationIndex(
        rep(valid_index$expected_mz, 2), valid_index, ppm = 0,
        ms1Spectrum = spectrum
    ))
    actual <- result[!unavailable, names(valid), drop = FALSE]
    rownames(actual) <- NULL
    expect_equal(actual, valid)
    expect_true(all(valid$formula_status == "parsed"))
    expect_warning(queryMZAnnotationIndex(observed, index, ppm = 0),
        "3 distinct matched reference formula"
    )
})

test_that("formula parsing is cached within a query without hiding other warnings", {
    target <- data.frame(formula = "CH3OH", exactmass = 32.026214748)
    index <- buildMZAnnotationIndex(target,
        adducts = c("M+H", "M+Na"), inferStructure = "never"
    )
    parser <- MetaboCoreUtils::countElements
    calls <- 0L
    local_mocked_bindings(countElements = function(x) {
        calls <<- calls + 1L
        parser(x)
    }, .package = "MetaboCoreUtils")
    observed <- rep(index$expected_mz, 3)
    queryMZAnnotationIndex(observed, index, ppm = 0,
        ms1Spectrum = data.frame(mz = index$expected_mz, intensity = 100)
    )
    expect_equal(calls, 1L)
    queryMZAnnotationIndex(index$expected_mz, index, ppm = 0)
    expect_equal(calls, 2L)

    local_mocked_bindings(countElements = function(x) {
        warning("Unexpected parser condition")
        parser(x)
    }, .package = "MetaboCoreUtils")
    expect_warning(unavailable <- .formula_counts("CH3OH"),
        "Unexpected parser condition"
    )
    expect_null(unavailable)
    local_mocked_bindings(countElements = function(x) {
        stop("Unexpected parser failure")
    }, .package = "MetaboCoreUtils")
    expect_error(.formula_counts("CH3OH"), "Unexpected parser failure")
})

test_that("disabled formula checks are quiet and empty results retain the schema", {
    index <- buildMZAnnotationIndex(
        data.frame(formula = "C5H9NO4R2", exactmass = 100),
        adducts = "M+H", inferStructure = "never"
    )
    expect_silent(result <- queryMZAnnotationIndex(index$expected_mz, index,
        ppm = 0, checkFormulaMass = FALSE, checkIsotopes = FALSE
    ))
    expect_identical(result$formula_status, "not_checked")
    expect_silent(empty <- queryMZAnnotationIndex(1000, index, ppm = 0))
    expect_identical(empty$formula_status, character())
    expect_identical(names(result), names(empty))
})

test_that("annotateSM stores inspectable formula status while keeping assay values", {
    index <- buildMZAnnotationIndex(
        data.frame(formula = "C5H9NO4R2", exactmass = 100),
        adducts = "M+H", inferStructure = "never"
    )
    counts <- matrix(c(5, 10), nrow = 1,
        dimnames = list("candidate", c("p1", "p2"))
    )
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(mz = index$expected_mz),
        spatialCoords = cbind(x = 1:2, y = c(0, 0))
    )
    expect_warning(result <- annotateSM(object,
        index = index, rawMzColumn = "mz", ppm = 0,
        returnOnlyAnnotated = FALSE, verbose = FALSE
    ), "1 distinct matched reference formula")
    candidates <- annotationCandidates(result)
    expect_identical(candidates$FormulaStatus, "unavailable")
    expect_true(is.na(candidates$ChemicalScore))
    expect_identical(SummarizedExperiment::assay(result), counts)
    empty <- annotateTable(data.frame(row_id = 1L, mz = 1000),
        index = index, ppm = 0, verbose = FALSE
    )
    expect_identical(empty$FormulaStatus, character())
})

test_that("water loss requires complete formulas even with structure capacity", {
    target <- data.frame(
        name = "labelled ethanol", formula = "C2H6O[59Co]",
        exactmass = 46.041864812, iso_smiles = "CCO"
    )
    index <- buildMZAnnotationIndex(target,
        adducts = c("M+H", "M+H-H2O"), minStructureScore = 0
    )
    expect_warning(result <- queryMZAnnotationIndex(index$expected_mz, index,
        ppm = 0
    ), "1 distinct matched reference formula")
    expect_identical(result$adduct, "M+H")
    expect_identical(result$formula_status, "unavailable")
    loss <- index$expected_mz[index$rules$name[index$rule_idx] == "M+H-H2O"]
    expect_warning(empty <- queryMZAnnotationIndex(loss, index, ppm = 0),
        "1 distinct matched reference formula"
    )
    expect_equal(nrow(empty), 0L)
})

test_that("decoy queries consolidate the distinct unavailable formulas", {
    index <- buildMZAnnotationIndex(
        data.frame(formula = "C5H9NO4R2", exactmass = c(99, 100, 101)),
        adducts = "M+H", inferStructure = "never", collapseIsomers = FALSE
    )
    warnings <- list()
    result <- withCallingHandlers(
        estimateAnnotationFDR(index$expected_mz[2], index,
            ppm = 1, shiftsDa = c(-1, 1)
        ),
        spamtp_formula_warning = function(condition) {
            warnings[[length(warnings) + 1L]] <<- condition
            invokeRestart("muffleWarning")
        }
    )
    expect_length(warnings, 1L)
    expect_identical(warnings[[1]]$formulas, "C5H9NO4R2")
    expect_true(all(result$best_scores[, -1] > 0))
})

test_that("halogen evidence is applied only with complete formula counts", {
    index <- buildMZAnnotationIndex(
        data.frame(formula = c("NaCl", "C2H4R", "C2H4[59Co]"),
            exactmass = c(57.95862196, 100, 200)
        ),
        polarity = "negative", adducts = "M+Cl", inferStructure = "never"
    )
    spectrum <- data.frame(mz = index$expected_mz, intensity = 100)
    expect_warning(result <- queryMZAnnotationIndex(index$expected_mz, index,
        ppm = 0, ms1Spectrum = spectrum
    ), "2 distinct matched reference formula")
    valid <- result$formula == "NaCl"
    expect_equal(result$isotope_score[valid], 0.2)
    expect_true(all(is.na(result$isotope_score[!valid])))
    expect_true(all(is.finite(result$score)))

    rule <- adductRules("positive")[1, , drop = FALSE]
    expect_identical(.formula_mass_score(100.9, "C2H4[59Co]", rule), NA_real_)
})

test_that("formula status describes each row independently of query order", {
    index <- buildMZAnnotationIndex(
        data.frame(formula = "C2H6O", exactmass = 46.041864812,
            iso_smiles = "CCO"
        ),
        adducts = c("M+H", "M+H-H2O"), minStructureScore = 0
    )
    forward <- queryMZAnnotationIndex(index$expected_mz, index, ppm = 0,
        checkFormulaMass = FALSE, checkIsotopes = FALSE
    )
    reverse <- queryMZAnnotationIndex(rev(index$expected_mz), index, ppm = 0,
        checkFormulaMass = FALSE, checkIsotopes = FALSE
    )
    expect_identical(forward, reverse)
    expect_identical(forward$formula_status[forward$adduct == "M+H"],
        "not_checked"
    )
    expect_identical(forward$formula_status[forward$adduct == "M+H-H2O"],
        "parsed"
    )
})
