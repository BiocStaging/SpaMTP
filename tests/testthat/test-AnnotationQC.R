qc_db <- function() {
    data.frame(
        formula = c("C5H8O5", "C6H12O6", "C4H9NO2"),
        exactmass = c(148.037173366, 180.063388104, 103.063328537),
        isomers = c("A", "B", "C"),
        isomers_names = c("2-hydroxyglutarate", "hexose", "GABA"),
        ramp_ids = c("RAMP_C_A", "RAMP_C_B", "RAMP_C_C")
    )
}

test_that("decoy FDR counts real matches and rejects shifts inside tolerance", {
    index <- buildMZAnnotationIndex(qc_db(), adducts = c("M+H", "M+Na"))
    targets <- index$expected_mz[1:3]
    fdr <- estimateAnnotationFDR(targets, index,
        ppm = 5,
        cutoffs = c(0, 0.5)
    )
    expect_s3_class(fdr, "spamtp_annotation_fdr")
    expect_equal(fdr$curve$target_peaks[fdr$curve$cutoff == 0], 3)
    expect_equal(fdr$curve$decoy_peaks, c(0, 0))
    expect_equal(fdr$curve$fdr[1], 0)
    expect_identical(colnames(fdr$best_scores)[1:2], c("mz", "target"))
    expect_length(fdr$shifts_da, 4)
    expect_error(
        estimateAnnotationFDR(targets, index,
            ppm = 5,
            shifts_da = c(-0.001, 0.001)
        ),
        "twice the tolerance"
    )
    expect_output(print(fdr), "mass-shift decoys")
})

test_that("decoys find chance matches when shifted onto another compound", {
    index <- buildMZAnnotationIndex(qc_db(), adducts = "M+H")
    gaba <- min(index$expected_mz) # GABA [M+H]+
    fdr <- estimateAnnotationFDR(gaba - 0.004, index,
        ppm = 5,
        shifts_da = c(-0.004, 0.004), cutoffs = 0
    )
    expect_equal(fdr$curve$target_peaks, 0)
    expect_equal(fdr$curve$decoy_peaks, 0.5)
    expect_true(is.na(fdr$curve$fdr))
})

test_that("decoy shifts keep adduct spacings and score duplicated m/z", {
    index <- buildMZAnnotationIndex(qc_db(), adducts = c("M+H", "M+Na"))
    mh <- index$expected_mz[1]
    fdr <- estimateAnnotationFDR(c(mh, mh), index, ppm = 5, cutoffs = 0)
    expect_equal(fdr$best_scores$target[1], fdr$best_scores$target[2])
    expect_gt(fdr$best_scores$target[1], 0)
    # a constant shift moves both ions of a pair by the same mass
    shift <- fdr$shifts_da[1]
    pair <- c(149.044450, 171.026394)
    expect_equal(diff(pair + shift), diff(pair))
})

test_that("score cutoff is the smallest one meeting the target FDR", {
    fdr <- structure(
        list(curve = data.frame(
            cutoff = c(0, 0.1, 0.5), target_peaks = c(100, 60, 10),
            decoy_peaks = c(90, 20, 1), fdr = c(0.9, 1 / 3, 0.1)
        )),
        class = "spamtp_annotation_fdr"
    )
    expect_equal(annotationScoreCutoff(fdr, 0.5), 0.1)
    expect_equal(annotationScoreCutoff(fdr, 0.1), 0.5)
    expect_true(is.na(annotationScoreCutoff(fdr, 0.1, min_peaks = 20)))
    expect_true(is.na(annotationScoreCutoff(fdr, 0.01)))
    expect_error(annotationScoreCutoff(list()), "estimateAnnotationFDR")
})

test_that("MALDI preset drops multiply charged and solvent adducts", {
    preset <- maldiAdductPreset("DHB")
    rules <- maldiMatrixRules("DHB", polarity = "positive")
    kept <- rules[rules$name %in% preset, ]
    expect_true(all(abs(kept$charge) == 1))
    expect_false(any(grepl("ACN|CH3OH|IsoProp|DMSO", preset)))
    expect_true(all(c("M+H", "M+Na", "M+K") %in% preset))
    expect_false("M+2H" %in% preset)
    expect_true(any(grepl("DHB", preset)))
    no_matrix <- maldiAdductPreset("DHB", include_matrix_adducts = FALSE)
    expect_false(any(grepl("DHB", no_matrix)))
})

test_that("matrix ion table reproduces known DHB ions", {
    ions <- matrixIonTable("DHB")
    mz_of <- function(label) ions$mz[ions$ion == label]
    expect_equal(mz_of("[M+H]+"), 155.033885, tolerance = 1e-6)
    expect_equal(mz_of("[M+H-H2O]+"), 137.023321, tolerance = 1e-6)
    expect_equal(mz_of("[2M+H-2H2O]+"), 273.039366, tolerance = 1e-6)
    expect_equal(mz_of("[M+Na]+"), 177.015830, tolerance = 1e-6)
    neg <- matrixIonTable("DHB", "negative")
    expect_true(all(grepl("\\]-$", neg$ion)))
    expect_equal(neg$mz[neg$ion == "[M-H]-"], 153.019332, tolerance = 1e-6)
    expect_error(matrixIonTable("fmp10"), "No matrix ion formula")
})

test_that("matrix peaks are flagged within tolerance only", {
    flags <- flagMatrixPeaks(c(155.0344, 273.0393, 200.1, 155.0400), "DHB",
        ppm = 5
    )
    expect_equal(flags$matrix_ion, c(TRUE, TRUE, FALSE, FALSE))
    expect_equal(flags$ion[1:2], c("[M+H]+", "[2M+H-2H2O]+"))
    expect_true(abs(flags$ppm_error[1]) < 5)
})

test_that("13C isotope peaks are flagged by spacing and intensity ratio", {
    mz <- c(399.1436, 400.1466, 410.2, 411.2033)
    # 400.1466 is the 13C peak of 399.1436; 411.2033 is spaced correctly but
    # is stronger than its putative parent, so it is not an isotope peak
    flags <- flagIsotopePeaks(mz, c(100, 14, 10, 50), ppm = 3)
    expect_equal(flags$isotope, c(FALSE, TRUE, FALSE, FALSE))
    expect_equal(flags$parent[2], 1L)
    expect_true(abs(flags$ppm_error[2]) < 3)
    expect_error(flagIsotopePeaks(mz, c(100, NA, 10, 50)), "finite")
    expect_error(flagIsotopePeaks(mz, c(100, -1, 10, 50)), "non-negative")
})

test_that("isotope flags require co-localisation when profiles are given", {
    set.seed(2)
    base <- abs(stats::rnorm(40)) + 1
    profiles <- rbind(
        base, 0.14 * base + stats::rnorm(40, sd = 0.01),
        base, abs(stats::rnorm(40)) * 0.14
    )
    mz <- c(399.1436, 400.1466, 500.2, 501.2033)
    flags <- flagIsotopePeaks(mz, profiles, ppm = 3, min_cor = 0.5)
    expect_equal(flags$isotope, c(FALSE, TRUE, FALSE, FALSE))
    expect_gt(flags$correlation[2], 0.9)
})

adduct_fixture <- function() {
    set.seed(3)
    m <- SpaMTP:::.spamtp_mass_constants()
    mh <- 149.044450
    list(
        m = m, mh = mh, mna = mh - m$proton + m$sodium_ion,
        base = abs(stats::rnorm(60)) + 1
    )
}

test_that("co-localised adduct peaks group, uncorrelated ones do not", {
    f <- adduct_fixture()
    profiles <- rbind(
        f$base, 0.7 * f$base + stats::rnorm(60, sd = 0.02),
        abs(stats::rnorm(60)), f$base * 1.0001
    )
    mz <- c(
        f$mh, f$mna, f$mh - f$m$proton + f$m$potassium_ion,
        f$mh * (1 + 3e-6)
    )
    g <- groupAdductPeaks(mz, profiles, ppm = 3, min_cor = 0.5, split_ppm = 5)
    expect_equal(g$groups$group[1], g$groups$group[2])
    expect_false(g$groups$group[3] == g$groups$group[1])
    expect_equal(g$groups$group[4], g$groups$group[1])
    expect_setequal(g$edges$relation, c("Na-H", "split"))
    # excluded peaks stay singletons
    gx <- groupAdductPeaks(mz, profiles,
        ppm = 3, min_cor = 0.5, split_ppm = 5,
        exclude = c(FALSE, TRUE, FALSE, FALSE)
    )
    expect_equal(gx$groups$group_size[2], 1L)
    expect_equal(gx$groups$group[1], gx$groups$group[4])
    expect_error(groupAdductPeaks(c(1, 2), c(1, 2)), "matrix")
    expect_error(groupAdductPeaks(mz, profiles, exclude = 9), "exclude")
})

test_that("group identities need support from two different adducts", {
    groups <- data.frame(
        peak = 1:3, mz = c(149.0444, 171.0264, 300.1),
        group = c(1, 1, 2), group_size = c(2, 2, 1)
    )
    cand <- data.frame(
        observed_mz = c(149.0444, 171.0264, 149.0444, 300.1),
        ramp_ids = c("RAMP_C_A", "RAMP_C_A", "RAMP_C_B; RAMP_C_C", "RAMP_C_D"),
        adduct = c("M+H", "M+Na", "M+K", "M+H"),
        score = c(0.4, 0.3, 0.6, 0.9)
    )
    a <- annotateAdductGroups(cand, groups)
    expect_equal(a$compounds$ramp_id, "RAMP_C_A")
    expect_equal(a$compounds$adducts, "M+H;M+Na")
    expect_equal(a$peaks$peak, 1:2)
    expect_equal(a$peaks$ramp_ids, c("RAMP_C_A", "RAMP_C_A"))
    expect_equal(a$peaks$n_identities, c(1L, 1L))
    expect_equal(nrow(annotateAdductGroups(cand, groups,
        min_score = 0.35
    )$peaks), 0)
    # candidate m/z off by rounding still matches its peak
    cand$observed_mz <- round(cand$observed_mz, 3)
    expect_equal(nrow(annotateAdductGroups(cand, groups,
        match_ppm = 5
    )$peaks), 2)
    cand$ramp_ids[1] <- NA
    expect_equal(nrow(annotateAdductGroups(cand, groups)$peaks), 0)
    expect_error(annotateAdductGroups(cand[, 1:2], groups), "candidates needs")
})

test_that("molecules chained into one group keep their own identities", {
    groups <- data.frame(
        peak = 1:4, mz = c(149.0444, 171.0264, 200.1, 222.08),
        group = 1, group_size = 4
    )
    cand <- data.frame(
        observed_mz = c(149.0444, 171.0264, 200.1, 222.08),
        ramp_ids = c("RAMP_C_A", "RAMP_C_A", "RAMP_C_E", "RAMP_C_E"),
        adduct = c("M+H", "M+Na", "M+H", "M+Na"), score = 0.5
    )
    a <- annotateAdductGroups(cand, groups)
    expect_setequal(a$compounds$ramp_id, c("RAMP_C_A", "RAMP_C_E"))
    expect_equal(
        a$peaks$ramp_ids,
        c("RAMP_C_A", "RAMP_C_A", "RAMP_C_E", "RAMP_C_E")
    )
})

test_that("database scope keeps compounds from the chosen pathway sources", {
    chem <- data.frame(ramp_id = c("RAMP_C_1", "RAMP_C_2", "RAMP_C_3"))
    ahp <- data.frame(
        rampId = c("RAMP_C_1", "RAMP_C_2"),
        pathwaySource = c("reactome", "pfocr")
    )
    expect_equal(annotationDatabaseScope(chem, ahp)$ramp_id, "RAMP_C_1")
    expect_equal(
        nrow(annotationDatabaseScope(chem, ahp,
            sources = c("reactome", "pfocr")
        )),
        2
    )
    expect_error(annotationDatabaseScope(data.frame(x = 1), ahp), "ramp_id")
})
