# Development assistance: see inst/COPYRIGHTS.

test_that("adduct ion masses reproduce MetaboCoreUtils::mass2mz()", {
    rules <- adductRules("both")
    mapped <- rules[!is.na(rules$metabocoreutils_name), , drop = FALSE]
    expect_gt(nrow(mapped), 0)
    neutral <- c(75.032028409, 148.037173366, 299.282428984, 885.549853)
    for (i in seq_len(nrow(mapped))) {
        adduct <- mapped$metabocoreutils_name[[i]]
        spamtp_mz <- (mapped$n_molecules[[i]] * neutral + mapped$mass_shift[[i]]) /
            abs(mapped$charge[[i]])
        expected <- as.numeric(MetaboCoreUtils::mass2mz(neutral, adduct = adduct))
        expect_equal(spamtp_mz, expected, tolerance = 1e-9, info = mapped$name[[i]])
    }

    mcu <- rbind(
        MetaboCoreUtils::adducts("positive"),
        MetaboCoreUtils::adducts("negative")
    )
    source <- mcu[match(mapped$metabocoreutils_name, mcu$name), , drop = FALSE]
    expect_equal(mapped$charge, as.integer(source$charge))
    expect_equal(mapped$polarity == "positive", source$positive)
})

test_that("only the documented rules are defined outside MetaboCoreUtils", {
    rules <- adductRules("both")
    local <- rules[is.na(rules$metabocoreutils_name), , drop = FALSE]
    expect_identical(local$name, c("M-H2O-H", "M"))

    water <- MetaboCoreUtils::calculateMass("H2O")[[1]]
    deprotonated <- rules$mass_shift[rules$name == "M-H"]
    expect_equal(local$mass_shift, c(deprotonated - water, 0), tolerance = 1e-12)
    expect_equal(local$charge, c(-1L, 1L))
    expect_equal(local$polarity, c("negative", "neutral"))
    expect_equal(local$loss_h, c(1L, 0L))
    expect_equal(local$added_charge, c(0L, 0L))

    matrix_rules <- maldiMatrixRules("DHB")
    matrix_only <- !matrix_rules$name %in% rules$name
    expect_true(any(matrix_only))
    expect_true(all(is.na(matrix_rules$metabocoreutils_name[matrix_only])))
})

test_that("polarity and complexity filtering keep SpaMTP rule names", {
    rules <- adductRules("both")
    expect_true(all(c("M+H", "M+Na", "2M+H", "M-H", "M+Cl", "M") %in% rules$name))
    for (mode in c("positive", "negative", "neutral")) {
        selected <- adductRules(mode)
        expect_true(all(selected$polarity == mode))
        expect_identical(selected$name, rules$name[rules$polarity == mode])
        expect_identical(names(selected), names(rules))
    }
    simple <- adductRules("negative", include_complex = FALSE)
    expect_identical(simple$name, c("M-H", "M+Cl", "M+Br"))
    expect_equal(simple$metabocoreutils_name, c("[M-H]-", "[M+Cl]-", "[M+Br]-"))
})
