#### Indexed and chemically constrained m/z annotation ########################

# Constants for SpaMTP-specific calculations (matrix products, QC spacings and
# isotope patterns). Standard adduct ion masses come from MetaboCoreUtils; see
# adductRules().
.spamtp_mass_constants <- function() {
    electron <- 0.000548579909065
    hydrogen <- 1.00782503223
    oxygen <- 15.99491461957

    list(
        electron = electron,
        proton = 1.007276466621,
        carbon = 12,
        hydrogen = hydrogen,
        nitrogen = 14.00307400443,
        oxygen = oxygen,
        sodium_ion = 22.9897692820 - electron,
        potassium_ion = 38.9637064864 - electron,
        chloride = 34.968852682 + electron,
        water = 2 * hydrogen + oxygen,
        carbon13_spacing = 1.00335483507,
        chlorine_m2_spacing = 1.99704991,
        bromine_m2_spacing = 1.9979535
    )
}

.make_adduct_rule <- function(
    name, polarity, nMolecules, charge,
    massShift, lossH = 0L,
    addedCharge = charge + lossH,
    complexity = "simple", baseAdduct = name,
    containsHalogen = FALSE,
    containsMetal = FALSE, prior = 1,
    waterLoss = 0L,
    metabocoreutilsName = NA_character_,
    carbonDelta = 0L, chlorineDelta = 0L, bromineDelta = 0L
) {
    data.frame(
        name = name,
        notation = paste0(
            "[", name, "]",
            ifelse(abs(charge) == 1, "", abs(charge)),
            ifelse(charge > 0, "+", "-")
        ),
        polarity = polarity,
        n_molecules = as.integer(nMolecules),
        charge = as.integer(charge),
        mass_shift = as.numeric(massShift),
        loss_h = as.integer(lossH),
        water_loss = as.integer(waterLoss),
        added_charge = as.integer(addedCharge),
        complexity = complexity,
        base_adduct = baseAdduct,
        contains_halogen = containsHalogen,
        contains_metal = containsMetal,
        prior = as.numeric(prior),
        metabocoreutils_name = as.character(metabocoreutilsName),
        carbon_delta = as.numeric(carbonDelta),
        chlorine_delta = as.numeric(chlorineDelta),
        bromine_delta = as.numeric(bromineDelta),
        stringsAsFactors = FALSE
    )
}

# SpaMTP short rule names and annotation metadata. `mcu` is the
# MetaboCoreUtils::adducts() name that supplies each ion mass.
.spamtp_adduct_map <- function() {
    utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
name        mcu            polarity loss_h water_loss complexity     base halogen metal prior
M+H         [M+H]+         positive 0      0          simple         M+H  FALSE   FALSE 1
M+H-H2O     [M+H-H2O]+     positive 0      1          neutral_loss   M+H  FALSE   FALSE 0.55
M+H-2H2O    [M+H-H4O2]+    positive 0      2          neutral_loss   M+H  FALSE   FALSE 0.35
M+NH4       [M+NH4]+       positive 0      0          salt           M+H  FALSE   FALSE 0.9
M+Na        [M+Na]+        positive 0      0          salt           M+Na FALSE   TRUE  0.9
M+K         [M+K]+         positive 0      0          salt           M+K  FALSE   TRUE  0.8
M+2H        [M+2H]2+       positive 0      0          multicharged   M+H  FALSE   FALSE 0.8
M+H+Na      [M+H+Na]2+     positive 0      0          multicharged   M+Na FALSE   TRUE  0.7
M+2Na       [M+2Na]2+      positive 0      0          multicharged   M+Na FALSE   TRUE  0.65
M+H+K       [M+H+K]2+      positive 0      0          multicharged   M+K  FALSE   TRUE  0.65
M+3H        [M+3H]3+       positive 0      0          multicharged   M+H  FALSE   FALSE 0.65
M+2H+Na     [M+2H+Na]3+    positive 0      0          multicharged   M+Na FALSE   TRUE  0.55
M+H+2Na     [M+H+Na2]3+    positive 0      0          multicharged   M+Na FALSE   TRUE  0.5
M+3Na       [M+Na3]3+      positive 0      0          multicharged   M+Na FALSE   TRUE  0.45
M+2Na-H     [M+2Na-H]+     positive 1      0          metal_exchange M+Na FALSE   TRUE  0.65
M+2K-H      [M+2K-H]+      positive 1      0          metal_exchange M+K  FALSE   TRUE  0.5
M+ACN+H     [M+C2H3N+H]+   positive 0      0          solvent        M+H  FALSE   FALSE 0.5
M+ACN+Na    [M+C2H3N+Na]+  positive 0      0          solvent        M+Na FALSE   TRUE  0.45
M+2ACN+H    [M+C4H6N2+H]+  positive 0      0          solvent        M+H  FALSE   FALSE 0.35
M+CH3OH+H   [M+CH4O+H]+    positive 0      0          solvent        M+H  FALSE   FALSE 0.45
M+IsoProp+H [M+C3H8O+H]+   positive 0      0          solvent        M+H  FALSE   FALSE 0.35
M+DMSO+H    [M+C2H6OS+H]+  positive 0      0          solvent        M+H  FALSE   FALSE 0.35
2M+H        [2M+H]+        positive 0      0          multimer       M+H  FALSE   FALSE 0.55
2M+NH4      [2M+NH4]+      positive 0      0          multimer       M+H  FALSE   FALSE 0.45
2M+Na       [2M+Na]+       positive 0      0          multimer       M+Na FALSE   TRUE  0.5
2M+K        [2M+K]+        positive 0      0          multimer       M+K  FALSE   TRUE  0.4
2M+ACN+H    [2M+C2H3N+H]+  positive 0      0          multimer       M+H  FALSE   FALSE 0.3
2M+ACN+Na   [2M+C2H3N+Na]+ positive 0      0          multimer       M+Na FALSE   TRUE  0.3
M-H         [M-H]-         negative 1      0          simple         M-H  FALSE   FALSE 1
M-2H        [M-2H]2-       negative 2      0          multicharged   M-H  FALSE   FALSE 0.8
M-3H        [M-3H]3-       negative 3      0          multicharged   M-H  FALSE   FALSE 0.6
M+Cl        [M+Cl]-        negative 0      0          salt           M-H  TRUE    FALSE 0.65
M+Br        [M+Br]-        negative 0      0          salt           M-H  TRUE    FALSE 0.55
M+Na-2H     [M+Na-2H]-     negative 2      0          metal_exchange M-H  FALSE   TRUE  0.6
M+K-2H      [M+K-2H]-      negative 2      0          metal_exchange M-H  FALSE   TRUE  0.5
M+FA-H      [M+CHO2]-      negative 0      0          solvent        M-H  FALSE   FALSE 0.55
M+Hac-H     [M+C2H3O2]-    negative 0      0          solvent        M-H  FALSE   FALSE 0.5
M+TFA-H     [M+C2F3O2]-    negative 0      0          solvent        M-H  TRUE    FALSE 0.35
M-H2O-H     NA             negative 1      0          neutral_loss   M-H  FALSE   FALSE 0.55
2M-H        [2M-H]-        negative 1      0          multimer       M-H  FALSE   FALSE 0.55
2M+FA-H     [2M+CHO2]-     negative 0      0          multimer       M-H  FALSE   FALSE 0.4
2M+Hac-H    [2M+C2H3O2]-   negative 0      0          multimer       M-H  FALSE   FALSE 0.35
3M-H        [3M-H]-        negative 1      0          multimer       M-H  FALSE   FALSE 0.3
M           NA             neutral  0      0          neutral        M    FALSE   FALSE 1
")
}

# Ion definitions for the two rules without a MetaboCoreUtils adduct: negative
# water loss (MetaboCoreUtils lists water losses for positive ions only) and the
# identity rule for neutral-mass input, which is not an ion.
.local_adduct_ions <- function(mcu) {
    data.frame(
        name = c("M-H2O-H", "M"),
        n_molecules = 1,
        charge = c(-1, 1),
        mass_shift = c(
            mcu$mass_add[mcu$name == "[M-H]-"] -
                MetaboCoreUtils::calculateMass("H2O")[[1]],
            0
        )
    )
}

#' Return the validated SpaMTP adduct rule table
#'
#' Ion masses are taken from [MetaboCoreUtils::adducts()]: each rule's
#' molecular multiplier and total numerator mass shift are the MetaboCoreUtils
#' `mass_multi` and `mass_add` scaled by the absolute charge, so the ion m/z
#' `(n_molecules * M + mass_shift) / abs(charge)` equals
#' [MetaboCoreUtils::mass2mz()]. SpaMTP keeps its short rule names (for
#' example `"M+H"`) and adds annotation priors and chemical constraints: the
#' number of protons removed from the analyte, the charge contributed by added
#' species, water losses, metal and halogen flags, a complexity class and a
#' base adduct. This allows rules to be checked before candidate generation.
#' Only negative-mode water loss (`"M-H2O-H"`, the MetaboCoreUtils `[M-H]-`
#' mass minus water) and the identity rule `"M"` for neutral-mass input are
#' defined locally.
#'
#' @param polarity One of `"both"`, `"positive"`, `"negative"`, or
#'   `"neutral"`.
#' @param includeComplex Include multimers, solvent clusters, metal-exchange,
#'   multiply charged ions, and neutral losses.
#'
#' @details Positive `[M+H-H2O]+` and `[M+H-2H2O]+` rules require
#'   structure-derived dehydration capacity. Alcohol OH groups and adjacent
#'   carbon-bound hydrogens must support distinct, valence-compatible losses.
#'   Their priors (0.55 and 0.35) are heuristics, not calibrated probabilities.
#'   Structural compatibility does not predict fragment abundance or establish
#'   identity; acquisition conditions and MS/MS evidence remain relevant.
#'
#' @return A data frame containing validated adduct rules. Column
#'   `metabocoreutils_name` names the MetaboCoreUtils adduct that supplies the
#'   ion mass (`NA` for the locally defined rules). `carbon_delta`,
#'   `chlorine_delta` and `bromine_delta` describe isotope-relevant atoms added
#'   or removed by the rule, in addition to `n_molecules` copies of the compound.
#' @examples
#' rules <- adductRules("positive")
#' rules[rules$name %in% c("M+H", "M+Na"), c("name", "mass_shift", "charge")]
#' @export
adductRules <- function(
    polarity = c("both", "positive", "negative", "neutral"),
    includeComplex = TRUE
) {
    polarity <- match.arg(polarity)
    map <- .spamtp_adduct_map()
    mcu <- rbind(
        MetaboCoreUtils::adducts("positive"),
        MetaboCoreUtils::adducts("negative")
    )

    ion <- mcu[match(map$mcu, mcu$name), ]
    z <- abs(ion$charge)
    ion <- data.frame(
        n_molecules = round(ion$mass_multi * z),
        charge = ion$charge, mass_shift = ion$mass_add * z
    )
    local <- .local_adduct_ions(mcu)
    ion[match(local$name, map$name), ] <- local[, names(ion)]
    missing <- is.na(ion$mass_shift)
    if (any(missing)) {
        stop(
            "MetaboCoreUtils does not define adduct(s): ",
            paste(map$mcu[missing], collapse = ", ")
        )
    }

    rules <- .make_adduct_rule(
        map$name, map$polarity, ion$n_molecules, ion$charge, ion$mass_shift,
        lossH = map$loss_h,
        addedCharge = ifelse(map$polarity == "neutral", 0L,
            ion$charge + map$loss_h
        ),
        complexity = map$complexity, baseAdduct = map$base,
        containsHalogen = map$halogen, containsMetal = map$metal,
        prior = map$prior, waterLoss = map$water_loss,
        metabocoreutilsName = map$mcu
    )
    # Include atoms contributed by solvents/adducts in isotope predictions.
    atom_counts <- function(x) {
        x[is.na(x) | !nzchar(x)] <- "C0"
        MetaboCoreUtils::countElements(x)
    }
    add <- atom_counts(mcu$formula_add)
    subtract <- atom_counts(mcu$formula_sub)
    isotope_elements <- c(C = "carbon_delta", Cl = "chlorine_delta", Br = "bromine_delta")
    matched <- match(map$mcu, mcu$name)
    for (symbol in names(isotope_elements)) {
        delta <- vapply(seq_along(add), function(i) {
            sum(add[[i]][names(add[[i]]) == symbol]) - sum(subtract[[i]][names(subtract[[i]]) == symbol])
        }, numeric(1))
        rules[[isotope_elements[[symbol]]]] <- ifelse(is.na(matched), 0, delta[matched])
    }

    rownames(rules) <- NULL
    if (polarity != "both") {
        rules <- rules[rules$polarity == polarity, , drop = FALSE]
    }
    if (!isTRUE(includeComplex)) {
        rules <- rules[rules$complexity %in% c("simple", "salt", "neutral"), , drop = FALSE]
    }
    rules
}

#' List MALDI matrix and on-tissue derivatization profiles
#'
#' The registry separates conventional MALDI matrices, reactive matrices, and
#' derivatization reagents. A registered profile does not necessarily imply
#' that a universal mass-shift rule is available: some reactions require a
#' coupling reagent, produce several products, or remain strongly dependent on
#' the analyte structure and acquisition conditions.
#'
#' @param matrix Optional matrix/reagent name or alias. When `NULL`, return the
#'   complete registry.
#'
#' @return A data frame describing supported matrix profiles and their current
#'   automatic-rule status.
#' @examples
#' maldiMatrixProfiles("DHB")
#' @export
maldiMatrixProfiles <- function(matrix = NULL) {
    profiles <- data.frame(
        matrix = c(
            "none", "dhb", "chca", "9aa", "dan", "norharmane",
            "fmp10", "fmp8", "fmp9", "dpp_tfb", "tmp_tfb", "n_mepyba",
            "dnph", "coniferyl_aldehyde", "dhba", "dhap", "girard_t",
            "girard_p", "2_picolylamine", "tmpa", "ampp_hatu", "tahs"
        ),
        display_name = c(
            "No specified MALDI matrix", "2,5-Dihydroxybenzoic acid (DHB)",
            "alpha-Cyano-4-hydroxycinnamic acid (CHCA)",
            "9-Aminoacridine (9-AA)", "1,5-Diaminonaphthalene (DAN)",
            "Norharmane", "FMP-10", "FMP-8", "FMP-9",
            "2,4-Diphenylpyrylium tetrafluoroborate (DPP-TFB)",
            "2,4,6-Trimethylpyrylium tetrafluoroborate (TMP-TFB)",
            "N-Methylpyridinium boronic acid (N-MePyBA)",
            "2,4-Dinitrophenylhydrazine (DNPH)", "Coniferyl aldehyde (CA)",
            "2,4-Dihydroxybenzaldehyde (DHBA)",
            "2,5-Dihydroxyacetophenone (DHAP)", "Girard reagent T",
            "Girard reagent P", "2-Picolylamine (2-PA)",
            "N,N,N-Trimethyl-2-(piperazin-1-yl)ethanaminium (TMPA)",
            "AMPP/HATU", "TAHS"
        ),
        category = c(
            "unspecified", rep("conventional_matrix", 5),
            rep("reactive_matrix", 10), rep("otcd_reagent", 6)
        ),
        default_polarity = c(
            "positive", "positive", "positive", "negative", "negative", "both",
            "positive", "positive", "positive", "positive", "positive",
            "positive", "positive", "positive", "positive", "positive",
            "positive", "positive", "positive", "positive", "positive", "positive"
        ),
        target_groups = c(
            "none", "none", "none", "none", "none", "none",
            "primary/secondary amine; phenolic hydroxyl",
            "primary/secondary amine; phenolic hydroxyl",
            "primary/secondary amine; phenolic hydroxyl", "primary amine",
            "primary amine", "catechol/1,2-diol", "aldehyde/ketone",
            "primary amine", "primary amine", "primary amine",
            "aldehyde/ketone", "aldehyde/ketone", "carboxylic acid",
            "carboxylic acid", "carboxylic acid; aldehyde", "catecholamine"
        ),
        automatic_rules = c(
            "standard", "validated_matrix_adduct", "validated_matrix_adduct",
            "standard", "standard", "standard", "validated_reactive_product",
            rep("profile_only", 15)
        ),
        reference = c(
            NA_character_, "10.1021/acs.analchem.0c04720",
            "10.1021/acs.analchem.0c04720", "10.1002/rcm.750",
            NA_character_, NA_character_, "10.1038/s41592-019-0551-3",
            "10.1038/s41592-019-0551-3", "10.1038/s41592-019-0551-3",
            "10.1007/s13361-015-1119-9", NA_character_,
            "10.1021/acs.analchem.8b03746", NA_character_, NA_character_, NA_character_,
            NA_character_, NA_character_, NA_character_,
            "10.1021/acs.analchem.6b01021", "10.1021/acs.analchem.0c02303",
            "10.1021/jasms.2c00336", NA_character_
        ),
        stringsAsFactors = FALSE
    )

    if (is.null(matrix)) {
        return(profiles)
    }
    canonical <- vapply(matrix, .normalise_maldi_matrix, character(1))
    profiles[match(canonical, profiles$matrix), , drop = FALSE]
}

.normalise_maldi_matrix <- function(matrix) {
    if (is.null(matrix) || !length(matrix)) {
        return("none")
    }
    if (length(matrix) != 1L || is.na(matrix) || !nzchar(trimws(matrix))) {
        stop("maldiMatrix must be one non-empty matrix/reagent name.")
    }
    key <- tolower(gsub("[^[:alnum:]]", "", matrix))
    aliases <- c(
        none = "none", unspecified = "none", nomatrix = "none",
        dhb = "dhb", `25dhb` = "dhb", dihydroxybenzoicacid = "dhb",
        chca = "chca", hcca = "chca", alphacyano4hydroxycinnamicacid = "chca",
        `9aa` = "9aa", `9aminoacridine` = "9aa",
        dan = "dan", `15dan` = "dan", diaminonaphthalene = "dan",
        norharmane = "norharmane", fmp10 = "fmp10", fmp8 = "fmp8",
        fmp9 = "fmp9", dpptfb = "dpp_tfb", dpp = "dpp_tfb",
        tmptfb = "tmp_tfb", tmp = "tmp_tfb", nmepyba = "n_mepyba",
        dnph = "dnph", coniferylaldehyde = "coniferyl_aldehyde",
        ca = "coniferyl_aldehyde", dhba = "dhba", dhap = "dhap",
        girardt = "girard_t", gt = "girard_t", girardp = "girard_p",
        gp = "girard_p", `2picolylamine` = "2_picolylamine",
        `2pa` = "2_picolylamine", tmpa = "tmpa", ampphatu = "ampp_hatu",
        ampp = "ampp_hatu", tahs = "tahs"
    )
    if (!key %in% names(aliases)) {
        stop(
            "Unknown MALDI matrix/reagent '", matrix, "'. See maldiMatrixProfiles()."
        )
    }
    unname(aliases[[key]])
}

.resolve_maldi_polarity <- function(polarity = NULL, maldiMatrix = NULL) {
    if (!is.null(polarity)) {
        return(match.arg(polarity, c("positive", "negative", "neutral")))
    }
    if (is.null(maldiMatrix)) {
        return("positive")
    }
    profile <- maldiMatrixProfiles(maldiMatrix)
    selected <- profile$default_polarity[[1]]
    if (identical(selected, "both")) {
        stop(
            profile$display_name[[1]],
            " supports both ion modes; supply polarity explicitly."
        )
    }
    selected
}

.decorate_matrix_rules <- function(
    rules, matrix, ruleClass, source,
    reactiveGroup = "none",
    minReactiveSites = 0L,
    siteCountColumn = NA_character_
) {
    n <- nrow(rules)
    rules$maldi_matrix <- rep_len(matrix, n)
    rules$rule_class <- rep_len(ruleClass, n)
    rules$rule_source <- rep_len(source, n)
    rules$reactive_group <- rep_len(reactiveGroup, n)
    rules$min_reactive_sites <- rep_len(as.integer(minReactiveSites), n)
    rules$site_count_column <- rep_len(siteCountColumn, n)
    rules
}

.matrix_specific_rules <- function(matrix, polarity) {
    m <- .spamtp_mass_constants()
    empty <- .decorate_matrix_rules(
        adductRules(polarity)[0, , drop = FALSE], matrix,
        "matrix_profile", NA_character_
    )
    if (polarity != "positive") {
        return(empty)
    }

    if (matrix == "fmp10") {
        # Net neutral atom additions; every observed product carries one
        # positive charge. For b, remove one H atom and one electron.
        retained_tag <- 20 * m$carbon + 14 * m$hydrogen + m$nitrogen
        rules <- rbind(
            .make_adduct_rule(
                "M+FMP10", "positive", 1, 1, retained_tag - m$electron,
                complexity = "matrix_derivative", prior = 0.8, carbonDelta = 20L
            ),
            .make_adduct_rule(
                "M+2FMP10a", "positive", 1, 1,
                2 * retained_tag - (m$carbon + 3 * m$hydrogen) - m$electron,
                complexity = "matrix_derivative", prior = 0.75, carbonDelta = 39L
            ),
            .make_adduct_rule(
                "M+2FMP10b", "positive", 1, 1,
                2 * retained_tag - m$hydrogen - m$electron,
                complexity = "matrix_derivative", prior = 0.7, carbonDelta = 40L
            )
        )
        rules <- .decorate_matrix_rules(
            rules, matrix, "reactive_product",
            "10.1038/s41592-019-0551-3; 10.1002/cmtd.202500062",
            "FMP-10-reactive site", c(1L, 2L, 2L), "fmp10_reactive_sites"
        )
        return(rules)
    }

    if (matrix == "dhb") {
        dhb <- 7 * m$carbon + 6 * m$hydrogen + 4 * m$oxygen
        dehydrated_dhb <- dhb - m$water
        rules <- rbind(
            .make_adduct_rule(
                "M+(DHB-H2O)+H", "positive", 1, 1,
                dehydrated_dhb + m$proton,
                complexity = "matrix_adduct",
                baseAdduct = "M+H", prior = 0.4, carbonDelta = 7L
            ),
            .make_adduct_rule(
                "M+(DHB-H2O)+Na", "positive", 1, 1,
                dehydrated_dhb + m$sodium_ion,
                complexity = "matrix_adduct",
                baseAdduct = "M+Na", containsMetal = TRUE, prior = 0.3, carbonDelta = 7L
            ),
            .make_adduct_rule(
                "M+2(DHB-H2O)+H", "positive", 1, 1,
                2 * dehydrated_dhb + m$proton,
                complexity = "matrix_adduct",
                baseAdduct = "M+H", prior = 0.2, carbonDelta = 14L
            )
        )
        return(.decorate_matrix_rules(
            rules, matrix, "matrix_adduct", "10.1021/acs.analchem.0c04720",
            "amine-enriched; parent-ion support required"
        ))
    }

    if (matrix == "chca") {
        chca <- 10 * m$carbon + 7 * m$hydrogen + m$nitrogen + 3 * m$oxygen
        rule <- .make_adduct_rule(
            "M+CHCA+Na", "positive", 1, 1, chca + m$sodium_ion,
            complexity = "matrix_adduct", baseAdduct = "M+Na",
            containsMetal = TRUE, prior = 0.25, carbonDelta = 10L
        )
        return(.decorate_matrix_rules(
            rule, matrix, "matrix_adduct", "10.1021/acs.analchem.0c04720",
            "amine-enriched; parent-ion support required"
        ))
    }
    empty
}

#' Select annotation rules from the MALDI matrix/reagent profile
#'
#' Conventional matrices select ordinary ion rules and, where validated,
#' matrix-metabolite adducts. Reactive matrices additionally select covalent
#' product rules. Profiles without a verified universal net transformation are
#' retained in the registry but do not generate speculative mass shifts.
#'
#' @param maldiMatrix Matrix or derivatization reagent name; aliases such as
#'   `"2,5-DHB"`, `"HCCA"`, `"9-AA"`, and `"FMP-10"` are accepted.
#' @param polarity `"positive"`, `"negative"`, or `"neutral"`. If `NULL`, use
#'   the profile default; profiles supporting both modes require an explicit
#'   choice.
#' @param includeStandard Include the validated ordinary adduct rules.
#' @param includeMatrixProducts Include validated matrix-adduct or reactive
#'   product rules.
#'
#' @return A validated adduct/reaction rule data frame accepted by
#'   [buildMZAnnotationIndex()].
#' @examples
#' head(maldiMatrixRules("DHB", polarity = "positive"))
#' @export
maldiMatrixRules <- function(
    maldiMatrix, polarity = NULL,
    includeStandard = TRUE,
    includeMatrixProducts = TRUE
) {
    matrix <- .normalise_maldi_matrix(maldiMatrix)
    profile <- maldiMatrixProfiles(matrix)
    polarity <- .resolve_maldi_polarity(polarity, matrix)

    standard <- .decorate_matrix_rules(
        adductRules(polarity)[if (isTRUE(includeStandard)) TRUE else FALSE, , drop = FALSE],
        matrix, "standard_adduct", "SpaMTP validated general adduct rules"
    )
    specific <- if (isTRUE(includeMatrixProducts)) {
        .matrix_specific_rules(matrix, polarity)
    } else {
        standard[0, , drop = FALSE]
    }

    if (isTRUE(includeMatrixProducts) &&
        profile$automatic_rules[[1]] == "profile_only") {
        warning(
            "The ", profile$display_name[[1]],
            " profile is registered, but SpaMTP does not yet bundle a universal ",
            "validated net product-mass rule. Standard adduct rules were selected; ",
            "supply a study-specific custom rule table when appropriate.",
            call. = FALSE
        )
    }
    .validate_adduct_rules(rbind(standard, specific))
}

.normalise_adduct_name <- function(x) {
    x <- gsub("\\s+", "", as.character(x))
    x <- sub("^\\[", "", x)
    x <- sub("\\](?:[0-9]*[+-]|[+-][0-9]*)$", "", x)
    aliases <- c("M+2K+H" = "M+2K-H")
    hit <- x %in% names(aliases)
    x[hit] <- unname(aliases[x[hit]])
    x
}

.validate_adduct_rules <- function(rules) {
    canonical_loss <- match(
        .normalise_adduct_name(rules$name),
        c("M+H-H2O", "M+H-2H2O")
    )
    canonical_loss[is.na(canonical_loss)] <- 0L
    if (!"water_loss" %in% names(rules)) rules$water_loss <- canonical_loss
    required <- c(
        "name", "polarity", "n_molecules", "charge", "mass_shift",
        "loss_h", "added_charge", "complexity", "base_adduct",
        "contains_halogen", "contains_metal", "prior"
    )
    missing_columns <- setdiff(required, names(rules))
    if (length(missing_columns)) {
        stop(
            "Adduct rule table is missing: ",
            paste(missing_columns, collapse = ", ")
        )
    }
    optional_defaults <- list(
        water_loss = 0L,
        carbon_delta = 0L, chlorine_delta = 0L, bromine_delta = 0L,
        maldi_matrix = "unspecified",
        rule_class = "custom_adduct",
        rule_source = "user-supplied rule",
        reactive_group = "none",
        min_reactive_sites = 0L,
        site_count_column = NA_character_
    )
    for (column in names(optional_defaults)) {
        if (!column %in% names(rules)) {
            rules[[column]] <- rep_len(optional_defaults[[column]], nrow(rules))
        }
    }
    if (any(!is.finite(rules$n_molecules) | rules$n_molecules < 1)) {
        stop("Every adduct rule must have nMolecules >= 1.")
    }
    if (any(!is.finite(rules$charge) | rules$charge == 0)) {
        stop("Every adduct rule must use a non-zero observed charge.")
    }
    non_neutral <- rules$polarity != "neutral"
    balanced <- rules$added_charge - rules$loss_h == rules$charge
    if (any(non_neutral & !balanced)) {
        bad <- rules$name[non_neutral & !balanced]
        stop(
            "Unbalanced adduct rule(s): ", paste(bad, collapse = ", "),
            ". addedCharge - lossH must equal charge."
        )
    }
    if (any(!is.finite(rules$mass_shift)) || any(!is.finite(rules$prior))) {
        stop("Adduct mass shifts and priors must be finite.")
    }
    if (any(!is.finite(rules$water_loss) | !rules$water_loss %in% 0:2)) {
        stop("waterLoss must contain integers from zero to two.")
    }
    if (any(canonical_loss > 0L & rules$water_loss != canonical_loss)) {
        stop("Canonical water-loss names must agree with water_loss.")
    }
    water <- rules$water_loss > 0L
    m <- .spamtp_mass_constants()
    # 1e-5 Da accepts MetaboCoreUtils' tabulated masses (rounded to ~1e-6 Da).
    if (any(water & (rules$polarity != "positive" |
        rules$n_molecules != 1L | rules$charge != 1L | rules$loss_h != 0L |
        abs(rules$mass_shift - (m$proton - rules$water_loss * m$water)) > 1e-5))) {
        stop("Water-loss rules require a singly protonated monomer and neutral water masses.")
    }
    rules$water_loss <- as.integer(rules$water_loss)
    rules$rule_class[water] <- "structure_constrained_neutral_loss"
    rules$reactive_group[water] <- "alcohol_OH_and_adjacent_carbon_H"
    rules$min_reactive_sites <- suppressWarnings(
        as.integer(rules$min_reactive_sites)
    )
    if (any(is.na(rules$min_reactive_sites) | rules$min_reactive_sites < 0L)) {
        stop("minReactiveSites must contain non-negative integers.")
    }
    rules
}

.first_existing_column <- function(x, choices, required = TRUE) {
    hit <- choices[choices %in% names(x)]
    if (length(hit)) {
        return(hit[[1]])
    }
    if (required) {
        stop(
            "Database must contain one of these columns: ",
            paste(choices, collapse = ", ")
        )
    }
    NULL
}

.paste_unique <- function(x) {
    x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
    if (length(x)) paste(x, collapse = "; ") else NA_character_
}

.normalise_metabolite_db <- function(
    db, collapseIsomers = TRUE,
    inferStructure = c("auto", "never", "always"),
    structureBackend = c("auto", "native"),
    BPPARAM = BiocParallel::SerialParam(),
    retainStructureRecords = FALSE
) {
    if (!is.data.frame(db)) {
        stop("db must be a data.frame or data.table.")
    }
    database_metadata <- attr(db, "spamtp_database", exact = TRUE)
    inferStructure <- match.arg(inferStructure)
    structureBackend <- match.arg(structureBackend)
    db <- as.data.frame(db, stringsAsFactors = FALSE)
    supplied_proton_col <- .first_existing_column(
        db,
        c("max_exchangeable_protons", "max_protons", "exchangeable_protons"),
        required = FALSE
    )
    supplied_proton_values <- if (is.null(supplied_proton_col)) {
        rep(NA_real_, nrow(db))
    } else {
        suppressWarnings(as.numeric(db[[supplied_proton_col]]))
    }
    supplied_source <- if ("proton_bound_source" %in% names(db)) {
        as.character(db$proton_bound_source)
    } else {
        rep(NA_character_, nrow(db))
    }
    supplied_proton_bound <- if (is.null(supplied_proton_col)) {
        rep(FALSE, nrow(db))
    } else {
        is.finite(supplied_proton_values) &
            (is.na(supplied_source) | supplied_source != "SMILES-inferred")
    }
    smiles_col <- .first_existing_column(
        db, c("iso_smiles", "canonical_smiles", "smiles", "SMILES"),
        required = FALSE
    )
    needs_precomputed <- !"structure_valid" %in% names(db) ||
        !any(!is.na(db$structure_valid) & db$structure_valid)
    if (inferStructure == "auto" && !is.null(smiles_col) && isTRUE(needs_precomputed) &&
        is.list(database_metadata) &&
        identical(database_metadata$source, "spamtpdb")) {
        precomputed <- tryCatch(
            .spamtp_db_resource(
                "smiles_features",
                version = database_metadata$version %||% "latest",
                source = "spamtpdb",
                localDir = database_metadata$local_dir,
                offline = isTRUE(database_metadata$offline)
            ),
            error = function(e) NULL
        )
        if (is.data.frame(precomputed)) {
            feature_key <- .first_existing_column(
                precomputed, c("smiles", "structure_smiles", "iso_smiles"),
                required = FALSE
            )
            if (!is.null(feature_key)) {
                feature_row <- match(
                    as.character(db[[smiles_col]]),
                    as.character(precomputed[[feature_key]])
                )
                feature_columns <- setdiff(names(precomputed), feature_key)
                for (column in feature_columns) {
                    if (!column %in% names(db)) {
                        db[[column]] <- precomputed[[column]][feature_row]
                    }
                }
            }
        }
    }
    if (!is.null(smiles_col) && inferStructure != "never") {
        has_smiles <- !is.na(db[[smiles_col]]) & nzchar(as.character(db[[smiles_col]]))
        should_infer <- inferStructure == "always" ||
            !"structure_valid" %in% names(db) ||
            any(has_smiles & is.na(db$structure_valid))
        if (isTRUE(should_infer)) {
            structure_count <- length(unique(as.character(
                db[[smiles_col]][!is.na(db[[smiles_col]]) & nzchar(db[[smiles_col]])]
            )))
            runtime_limit <- getOption("SpaMTP.max_runtime_smiles", 5000L)
            runtime_limit <- suppressWarnings(as.integer(runtime_limit))
            if (length(runtime_limit) != 1L || is.na(runtime_limit) ||
                runtime_limit < 0L) {
                runtime_limit <- 5000L
            }
            if (inferStructure == "auto" &&
                structure_count > runtime_limit) {
                warning(
                    "Skipping runtime SMILES decomposition for ",
                    format(structure_count, big.mark = ","), " unique structures ",
                    "(SpaMTP.max_runtime_smiles = ", runtime_limit, "). Use a ",
                    "SpaMTPdb resource with precomputed structural fields, set ",
                    "inferStructure = 'always', or raise the option explicitly.",
                    call. = FALSE
                )
            } else {
                db <- annotateSMILESStructure(
                    db,
                    smilesColumn = smiles_col, backend = structureBackend,
                    overwrite = inferStructure == "always", strict = FALSE,
                    BPPARAM = BPPARAM
                )
            }
        }
    }
    formula_col <- .first_existing_column(db, c("formula", "mol_formula"))
    mass_col <- .first_existing_column(
        db, c(
            "exactmass",
            "monoisotop_mass",
            "monoisotopic_mass",
            "neutral_mass"
        )
    )
    id_col <- .first_existing_column(
        db, c("isomers", "chem_source_id", "compound_id", "id", "ramp_id"),
        required = FALSE
    )
    name_col <- .first_existing_column(
        db, c("isomers_names", "common_name", "name"),
        required = FALSE
    )
    inchikey_col <- .first_existing_column(
        db, c("isomers_inchikey", "inchi_key", "inchikey"),
        required = FALSE
    )
    ramp_col <- .first_existing_column(
        db,
        c("ramp_id", "rampId"),
        required = FALSE
    )
    proton_col <- .first_existing_column(
        db,
        c("max_exchangeable_protons", "max_protons", "exchangeable_protons"),
        required = FALSE
    )
    site_columns <- c(
        "fmp10_reactive_sites", "primary_amine_sites", "secondary_amine_sites",
        "phenolic_hydroxyl_sites", "catechol_sites", "carbonyl_sites",
        "carboxyl_sites", "hydroxyl_sites", "alcohol_hydroxyl_sites",
        "ketone_sites", "aldehyde_sites", "tertiary_amine_sites",
        "aromatic_nitrogen_sites", "amide_sites", "phosphate_acid_sites",
        "sulfonic_acid_sites", "thiol_sites", "acidic_sites",
        "weak_acidic_sites", "basic_sites", "proton_acceptor_sites",
        "alkali_binding_sites", "alkali_chelation_motifs",
        "positive_mode_score", "negative_mode_score",
        "alkali_affinity_score", "neutral_mass_score",
        "dehydration_sites", "max_water_losses"
    )

    n <- nrow(db)
    valid_structure <- if ("structure_valid" %in% names(db)) {
        !is.na(db$structure_valid) & as.logical(db$structure_valid)
    } else {
        rep(FALSE, n)
    }
    value_or_na <- function(column) {
        if (is.null(column)) rep(NA_character_, n) else as.character(db[[column]])
    }
    compounds <- data.frame(
        formula = as.character(db[[formula_col]]),
        exactmass = suppressWarnings(as.numeric(db[[mass_col]])),
        isomers = value_or_na(id_col),
        isomers_inchikey = value_or_na(inchikey_col),
        isomers_names = value_or_na(name_col),
        ramp_ids = value_or_na(ramp_col),
        smiles = if (is.null(smiles_col)) {
            rep(NA_character_, n)
        } else {
            as.character(db[[smiles_col]])
        },
        structure_valid = if ("structure_valid" %in% names(db)) {
            as.logical(db$structure_valid)
        } else {
            rep(NA, n)
        },
        structure_backend = if ("structure_backend" %in% names(db)) {
            as.character(db$structure_backend)
        } else {
            rep(NA_character_, n)
        },
        structure_evidence = if ("structure_evidence" %in% names(db)) {
            as.character(db$structure_evidence)
        } else {
            rep(NA_character_, n)
        },
        positive_site_atoms = if ("positive_site_atoms" %in% names(db)) {
            as.character(db$positive_site_atoms)
        } else {
            rep(NA_character_, n)
        },
        negative_site_atoms = if ("negative_site_atoms" %in% names(db)) {
            as.character(db$negative_site_atoms)
        } else {
            rep(NA_character_, n)
        },
        alkali_site_atoms = if ("alkali_site_atoms" %in% names(db)) {
            as.character(db$alkali_site_atoms)
        } else {
            rep(NA_character_, n)
        },
        proton_bound_source = ifelse(
            supplied_proton_bound,
            "provided",
            ifelse(valid_structure, "SMILES-inferred", "unavailable")
        ),
        max_exchangeable_protons = if (is.null(proton_col)) {
            rep(NA_real_, n)
        } else {
            inferred_or_primary <- suppressWarnings(as.numeric(db[[proton_col]]))
            inferred_or_primary[supplied_proton_bound] <-
                supplied_proton_values[supplied_proton_bound]
            inferred_or_primary
        },
        stringsAsFactors = FALSE
    )
    for (column in site_columns) {
        compounds[[column]] <- if (column %in% names(db)) {
            suppressWarnings(as.numeric(db[[column]]))
        } else {
            rep(NA_real_, n)
        }
    }
    keep <- is.finite(compounds$exactmass) & compounds$exactmass > 0 &
        !is.na(compounds$formula) & nzchar(compounds$formula)
    compounds <- compounds[keep, , drop = FALSE]
    structure_records <- if (retainStructureRecords) compounds else NULL
    structure_record_rows <- if (retainStructureRecords) {
        as.list(seq_len(nrow(compounds)))
    } else {
        NULL
    }

    if (isTRUE(collapseIsomers) && nrow(compounds)) {
        # Keep existing adduct groups stable; water losses resolve original records.
        grouping_columns <- setdiff(
            site_columns,
            c("dehydration_sites", "max_water_losses")
        )
        group_key <- paste(
            compounds$formula,
            format(
                compounds$exactmass,
                digits = 17,
                scientific = FALSE,
                trim = TRUE
            ),
            ifelse(is.na(compounds$max_exchangeable_protons), "NA",
                compounds$max_exchangeable_protons
            ),
            compounds$proton_bound_source,
            do.call(
                paste,
                c(
                    lapply(compounds[grouping_columns], function(x) {
                        ifelse(is.na(x), "NA", format(x, trim = TRUE))
                    }),
                    sep = ":"
                )
            ),
            sep = "\r"
        )
        if (anyDuplicated(group_key)) {
            grouped_rows <- split(seq_len(nrow(compounds)), group_key)
            if (retainStructureRecords) structure_record_rows <- unname(grouped_rows)
            first_row <- vapply(grouped_rows, function(i) i[[1]], integer(1))
            collapsed_site_values <- lapply(
                compounds[site_columns], function(x) x[first_row]
            )
            collapse_column <- function(column) {
                vapply(
                    grouped_rows,
                    function(i) .paste_unique(column[i]),
                    character(1)
                )
            }
            compounds <- data.frame(
                formula = compounds$formula[first_row],
                exactmass = compounds$exactmass[first_row],
                isomers = collapse_column(compounds$isomers),
                isomers_inchikey = collapse_column(compounds$isomers_inchikey),
                isomers_names = collapse_column(compounds$isomers_names),
                ramp_ids = collapse_column(compounds$ramp_ids),
                smiles = collapse_column(compounds$smiles),
                structure_valid = compounds$structure_valid[first_row],
                structure_backend = collapse_column(compounds$structure_backend),
                structure_evidence = collapse_column(compounds$structure_evidence),
                positive_site_atoms = collapse_column(compounds$positive_site_atoms),
                negative_site_atoms = collapse_column(compounds$negative_site_atoms),
                alkali_site_atoms = collapse_column(compounds$alkali_site_atoms),
                proton_bound_source = compounds$proton_bound_source[first_row],
                max_exchangeable_protons =
                    compounds$max_exchangeable_protons[first_row],
                stringsAsFactors = FALSE
            )
            for (column in site_columns) {
                compounds[[column]] <- collapsed_site_values[[column]]
            }
        }
    }
    rownames(compounds) <- NULL
    attr(compounds, "structure_records") <- structure_records
    attr(compounds, "structure_record_rows") <- structure_record_rows
    compounds
}

# Resolve only mass-matched original records, so an eligible isomer cannot
# confer dehydration eligibility on other members of a collapsed group.
.resolve_dehydration <- function(
    compounds, rules, inferStructure, cache, formulaCache = NULL
) {
    n <- nrow(compounds)
    for (column in c("dehydration_sites", "max_water_losses")) {
        if (!column %in% names(compounds)) compounds[[column]] <- rep(NA_real_, n)
    }
    status <- rep("not_required", n)
    keep <- rep(TRUE, n)
    for (j in which(rules$water_loss > 0L)) {
        smiles <- compounds$smiles[[j]]
        missing_prior <- !"positive_mode_score" %in% names(compounds) ||
            !is.finite(compounds$positive_mode_score[[j]])
        if (inferStructure != "never" &&
            (!isTRUE(compounds$structure_valid[[j]]) || missing_prior) &&
            !is.na(smiles) && nzchar(smiles)) {
            full_key <- paste0("full\n", smiles)
            if (!exists(full_key, envir = cache, inherits = FALSE)) {
                assign(full_key, .deconvolve_smiles_one(smiles), envir = cache)
            }
            features <- get(full_key, envir = cache, inherits = FALSE)
            for (column in intersect(names(features), names(compounds))) {
                if (is.na(compounds[[column]][[j]]) || column == "structure_valid") {
                    compounds[[column]][[j]] <- features[[column]]
                }
            }
            assign(smiles, features[c("dehydration_sites", "max_water_losses")],
                envir = cache
            )
        }
        capacity <- compounds$max_water_losses[[j]]
        if (!is.finite(capacity) &&
            isTRUE(compounds$structure_valid[[j]]) &&
            "alcohol_hydroxyl_sites" %in% names(compounds) &&
            isTRUE(compounds$alcohol_hydroxyl_sites[[j]] == 0)) {
            capacity <- 0L
            compounds$max_water_losses[[j]] <- 0L
            compounds$dehydration_sites[[j]] <- 0L
        }
        if (!is.finite(capacity) && inferStructure != "never") {
            if (!is.na(smiles) && nzchar(smiles)) {
                if (!exists(smiles, envir = cache, inherits = FALSE)) {
                    assign(
                        smiles,
                        .smiles_dehydration_capacity(smiles),
                        envir = cache
                    )
                }
                details <- get(smiles, envir = cache, inherits = FALSE)
                capacity <- details$max_water_losses
                compounds$max_water_losses[[j]] <- capacity
                compounds$dehydration_sites[[j]] <- details$dehydration_sites
            }
        }
        required <- rules$water_loss[[j]]
        atoms_available <- if ("formula" %in% names(compounds)) {
            isTRUE(.element_count(
                compounds$formula[[j]], "O", formulaCache
            ) >= required) && isTRUE(.element_count(
                compounds$formula[[j]], "H", formulaCache
            ) >= 2L * required)
        } else {
            TRUE
        }
        keep[[j]] <- is.finite(capacity) && capacity >= required && atoms_available
        status[[j]] <- if (!is.finite(capacity)) {
            "unavailable"
        } else if (keep[[j]]) {
            "structure-compatible"
        } else {
            "insufficient"
        }
        # A successfully checked motif supplies structure evidence even when a
        # large, unannotated library skipped full runtime decomposition.
        if (keep[[j]] && !isTRUE(compounds$structure_valid[[j]])) {
            compounds$structure_valid[[j]] <- TRUE
            compounds$structure_evidence[[j]] <-
                "Alcohol OH and adjacent carbon H support the requested water loss"
        }
    }
    list(compounds = compounds, keep = keep, status = status)
}

.adduct_structure_prior <- function(compounds, rule) {
    n <- nrow(compounds)
    score <- rep(1, n)
    status <- rep("unavailable", n)
    rationale <- rep("No parsed SMILES; structure prior not applied", n)
    parsed <- !is.na(compounds$structure_valid) & compounds$structure_valid
    if (!any(parsed)) {
        return(list(score = score, status = status, rationale = rationale))
    }

    value <- function(column, default = 0) {
        x <- if (column %in% names(compounds)) compounds[[column]] else rep(NA_real_, n)
        x[!is.finite(x)] <- default
        x
    }
    positive <- value("positive_mode_score", 0.25)
    negative <- value("negative_mode_score", 0.15)
    alkali <- value("alkali_affinity_score", 0.15)
    neutral <- value("neutral_mass_score", 1)
    acidic <- value("acidic_sites", 0)
    weak_acidic <- value("weak_acidic_sites", 0)
    acceptors <- value("proton_acceptor_sites", 0)
    basic <- value("basic_sites", 0)
    donors <- value("alkali_binding_sites", 0)

    selected <- if (rule$polarity[[1]] == "neutral") {
        neutral
    } else if (isTRUE(rule$contains_metal[[1]])) {
        alkali
    } else if (rule$polarity[[1]] == "positive") {
        positive
    } else if (rule$loss_h[[1]] > 0L) {
        negative
    } else {
        pmin(0.9, 0.35 + 0.08 * (acceptors + donors))
    }

    if (rule$loss_h[[1]] > 0L) {
        preferred_capacity <- acidic
        weak_capacity <- acidic + weak_acidic
        insufficient <- preferred_capacity < rule$loss_h[[1]]
        selected[insufficient & weak_capacity >= rule$loss_h[[1]]] <-
            selected[insufficient & weak_capacity >= rule$loss_h[[1]]] * 0.55
        selected[weak_capacity < rule$loss_h[[1]]] <-
            selected[weak_capacity < rule$loss_h[[1]]] * 0.05
    }
    if (rule$polarity[[1]] == "positive" && abs(rule$charge[[1]]) > 1L &&
        !isTRUE(rule$contains_metal[[1]])) {
        capacity <- basic + acceptors
        selected[capacity < abs(rule$charge[[1]])] <-
            selected[capacity < abs(rule$charge[[1]])] * 0.25
    }

    score[parsed] <- pmax(0.01, pmin(1, selected[parsed]))
    status[parsed] <- ifelse(score[parsed] >= 0.5, "structure-supported",
        "structure-weak"
    )
    family <- if (rule$polarity[[1]] == "neutral") {
        "neutral-mass state"
    } else if (isTRUE(rule$contains_metal[[1]])) {
        "alkali/metal coordination"
    } else if (rule$polarity[[1]] == "positive") {
        "positive-mode proton/cation acceptance"
    } else if (rule$loss_h[[1]] > 0L) {
        "negative-mode deprotonation"
    } else {
        "negative-mode anion attachment"
    }
    rationale[parsed] <- paste0(
        family, ": ", compounds$structure_evidence[parsed]
    )
    list(score = score, status = status, rationale = rationale)
}

.matrix_structure_compatibility <- function(compounds, matrix) {
    n <- nrow(compounds)
    matrix <- .normalise_maldi_matrix(matrix)
    profile <- maldiMatrixProfiles(matrix)
    category <- profile$category[[1]]
    target <- profile$target_groups[[1]]
    if (!category %in% c("reactive_matrix", "otcd_reagent")) {
        return(list(
            sites = rep(NA_real_, n), status = rep("not_reactive", n),
            target = rep(target, n)
        ))
    }
    value <- function(column) {
        x <- if (column %in% names(compounds)) compounds[[column]] else rep(NA_real_, n)
        suppressWarnings(as.numeric(x))
    }
    primary <- value("primary_amine_sites")
    secondary <- value("secondary_amine_sites")
    phenol <- value("phenolic_hydroxyl_sites")
    catechol <- value("catechol_sites")
    aldehyde <- value("aldehyde_sites")
    ketone <- value("ketone_sites")
    carboxyl <- value("carboxyl_sites")
    sites <- switch(matrix,
        fmp10 = primary + secondary + phenol,
        fmp8 = primary + secondary + phenol,
        fmp9 = primary + secondary + phenol,
        dpp_tfb = primary,
        tmp_tfb = primary,
        # Catechol is explicit. Non-adjacent alcohol counts must not be treated as
        # proof of the 1,2-diol motif targeted by N-MePyBA.
        n_mepyba = catechol,
        dnph = aldehyde + ketone,
        coniferyl_aldehyde = primary,
        dhba = primary,
        dhap = primary,
        girard_t = aldehyde + ketone,
        girard_p = aldehyde + ketone,
        `2_picolylamine` = carboxyl,
        tmpa = carboxyl,
        ampp_hatu = carboxyl + aldehyde,
        tahs = pmin(catechol, as.numeric(primary + secondary > 0)),
        rep(NA_real_, n)
    )
    parsed <- !is.na(compounds$structure_valid) & compounds$structure_valid
    status <- rep("unknown", n)
    status[parsed & is.finite(sites) & sites > 0] <- "compatible-target"
    status[parsed & is.finite(sites) & sites <= 0] <- "no-compatible-target"
    list(sites = sites, status = status, target = rep(target, n))
}

#' Predict a structure-aware adduct search space from SMILES
#'
#' This function applies the same functional-group and ion-mode logic used by
#' [buildMZAnnotationIndex()] before any observed m/z is supplied. It is useful
#' for auditing which protonation, deprotonation, alkali-binding, or reactive
#' matrix hypotheses SpaMTP will retain for a structure.
#'
#' @param smiles Character vector of SMILES strings.
#' @param polarity `"positive"`, `"negative"`, or `"neutral"`. When `NULL`,
#'   use the matrix-profile default or positive mode.
#' @param maldiMatrix Optional MALDI matrix/reagent profile.
#' @param rules Optional custom rule table. It cannot be combined with
#'   `maldiMatrix`.
#' @param minStructureScore Minimum rule-specific structural prior marked as
#'   retained.
#' @param backend,BPPARAM Passed to [deconvolveSMILES()].
#'
#' @return A ranked data frame with one row per structure/rule combination.
#' @examples
#' predicted <- predictAdductsFromSMILES("CC(O)C(O)=O", polarity = "negative")
#' head(predicted[, c("adduct", "structure_score", "retained")])
#' @export
predictAdductsFromSMILES <- function(
    smiles, polarity = NULL, maldiMatrix = NULL, rules = NULL,
    minStructureScore = 0.05,
    backend = c("auto", "native"),
    BPPARAM = BiocParallel::SerialParam()
) {
    backend <- match.arg(backend)
    polarity <- .resolve_maldi_polarity(polarity, maldiMatrix)
    if (is.null(rules)) {
        rules <- if (is.null(maldiMatrix)) {
            adductRules(polarity)
        } else {
            maldiMatrixRules(maldiMatrix, polarity = polarity)
        }
    } else {
        if (!is.null(maldiMatrix)) {
            stop("Supply either maldiMatrix or rules, not both.")
        }
        rules <- .validate_adduct_rules(as.data.frame(rules))
        rules <- rules[rules$polarity == polarity, , drop = FALSE]
    }
    rules <- .validate_adduct_rules(rules)
    if (!nrow(rules)) stop("No rules are available for the selected polarity.")
    if (length(minStructureScore) != 1L || !is.finite(minStructureScore) ||
        minStructureScore < 0 || minStructureScore > 1) {
        stop("minStructureScore must be a single number between zero and one.")
    }

    structures <- deconvolveSMILES(
        smiles,
        backend = backend, strict = FALSE, BPPARAM = BPPARAM
    )
    matrix_profile <- if (is.null(maldiMatrix)) {
        "none"
    } else {
        .normalise_maldi_matrix(maldiMatrix)
    }
    matrix_compatibility <- .matrix_structure_compatibility(
        structures, matrix_profile
    )
    parts <- vector("list", nrow(structures) * nrow(rules))
    k <- 0L
    for (i in seq_len(nrow(structures))) {
        compound <- structures[i, , drop = FALSE]
        for (j in seq_len(nrow(rules))) {
            k <- k + 1L
            rule <- rules[j, , drop = FALSE]
            prior <- .adduct_structure_prior(compound, rule)
            site_status <- "not_required"
            required <- rule$min_reactive_sites[[1]]
            site_column <- rule$site_count_column[[1]]
            if (required > 0L) {
                if (is.na(site_column) || !site_column %in% names(compound) ||
                    !is.finite(compound[[site_column]][[1]])) {
                    site_status <- "unknown"
                } else if (compound[[site_column]][[1]] >= required) {
                    site_status <- "verified"
                } else {
                    site_status <- "insufficient"
                }
            }
            site_score <- if (site_status == "unknown") 0.25 else 1
            proton_bound_status <- "not_required"
            if (rule$loss_h[[1]] > 0L && isTRUE(compound$structure_valid[[1]])) {
                preferred <- compound$acidic_sites[[1]]
                weak <- preferred + compound$weak_acidic_sites[[1]]
                proton_bound_status <- if (preferred >= rule$loss_h[[1]]) {
                    "preferred-sites"
                } else if (preferred > 0) {
                    "insufficient"
                } else if (weak >= rule$loss_h[[1]]) {
                    "weak-sites-only"
                } else {
                    "insufficient"
                }
            }
            retained <- prior$score[[1]] >= minStructureScore &&
                site_status != "insufficient" && proton_bound_status != "insufficient"
            dehydration <- .resolve_dehydration(
                compound, rule, "never", new.env(parent = emptyenv())
            )
            retained <- retained && dehydration$keep[[1]]
            parts[[k]] <- data.frame(
                smiles = compound$smiles,
                polarity = polarity,
                maldi_matrix = rule$maldi_matrix,
                adduct = rule$name,
                notation = rule$notation,
                rule_class = rule$rule_class,
                rule_prior = rule$prior,
                water_loss = rule$water_loss,
                dehydration_sites = compound$dehydration_sites,
                max_water_losses = compound$max_water_losses,
                dehydration_status = dehydration$status,
                structure_score = prior$score[[1]],
                reactive_site_status = site_status,
                proton_bound_status = proton_bound_status,
                matrix_target_groups = matrix_compatibility$target[[i]],
                matrix_target_sites = matrix_compatibility$sites[[i]],
                matrix_target_status = matrix_compatibility$status[[i]],
                combined_prior = rule$prior * prior$score[[1]] * site_score,
                retained = retained,
                structure_status = prior$status[[1]],
                structure_rule_evidence = prior$rationale[[1]],
                structure_evidence = compound$structure_evidence,
                stringsAsFactors = FALSE
            )
        }
    }
    result <- do.call(rbind, parts)
    result <- result[order(
        match(result$smiles, unique(as.character(smiles))),
        !result$retained, -result$combined_prior, result$adduct
    ), , drop = FALSE]
    rownames(result) <- NULL
    result
}

#' Build a reusable indexed metabolite annotation search space
#'
#' Expected adduct m/z values are generated once, sorted, and queried by binary
#' interval search. This is the one-dimensional equivalent of an interval tree
#' and avoids scanning every database row for every observed peak.
#'
#' @param db Metabolite database. Both the legacy SpaMTP five-column/wide
#'   database and the RaMP `chem_props` schema are supported.
#' @param polarity `"positive"`, `"negative"`, or `"neutral"`. When `NULL`,
#'   use the matrix-profile default, or positive mode when no profile is given.
#' @param adducts Optional character vector of adduct names or bracketed
#'   notations. When `NULL`, use the complete automatically selected rule
#'   space; this is a filter, not a compulsory input.
#' @param rules Optional custom rule table. See [adductRules()].
#' @param maldiMatrix Optional MALDI matrix or derivatization reagent profile.
#'   When supplied and `rules` is `NULL`, [maldiMatrixRules()] automatically
#'   selects standard and validated matrix-specific rules. `adducts` remains an
#'   optional filter on that selected rule space.
#' @param collapseIsomers Collapse records sharing formula, exact mass, and
#'   proton bound before indexing.
#' @param inferStructure `"auto"` joins precomputed SpaMTPdb features or
#'   derives them for at most `getOption("SpaMTP.max_runtime_smiles", 5000)`
#'   unique SMILES, `"never"` disables inference, and `"always"` forces
#'   runtime parsing and replaces precomputed structural fields. `"never"`
#'   and `"always"` do not fetch a precomputed structure table.
#' @param structureBackend SMILES parser used by [deconvolveSMILES()].
#' @param BPPARAM BiocParallel parameters for runtime SMILES parsing.
#' @param minStructureScore Minimum rule-specific structural prior retained
#'   before m/z indexing. Set to zero to rank without structure-based pruning.
#'
#' @return An object of class `spamtp_mz_index`.
#' @examples
#' data(spamtp_example)
#' index <- buildMZAnnotationIndex(spamtp_example$metabolites, adducts = "M+H")
#' queryMZAnnotationIndex(149.04445, index, ppm = 5)
#' @export
buildMZAnnotationIndex <- function(
    db, polarity = NULL,
    adducts = NULL, rules = NULL,
    maldiMatrix = NULL,
    collapseIsomers = TRUE,
    inferStructure = c("auto", "never", "always"),
    structureBackend = c("auto", "native"),
    BPPARAM = BiocParallel::SerialParam(),
    minStructureScore = 0.05
) {
    inferStructure <- match.arg(inferStructure)
    structureBackend <- match.arg(structureBackend)
    if (length(minStructureScore) != 1L || !is.finite(minStructureScore) ||
        minStructureScore < 0 || minStructureScore > 1) {
        stop("minStructureScore must be a single number between zero and one.")
    }
    polarity <- .resolve_maldi_polarity(polarity, maldiMatrix)
    if (is.null(rules)) {
        rules <- if (is.null(maldiMatrix)) {
            adductRules(polarity = polarity)
        } else {
            maldiMatrixRules(maldiMatrix, polarity = polarity)
        }
    } else {
        if (!is.null(maldiMatrix)) {
            stop(
                "Supply either maldiMatrix for automatic selection or an explicit ",
                "rules table, not both."
            )
        }
        rules <- as.data.frame(rules, stringsAsFactors = FALSE)
        rules <- rules[rules$polarity == polarity, , drop = FALSE]
    }
    rules <- .validate_adduct_rules(rules)
    family_rules <- rules

    if (!is.null(adducts)) {
        requested <- unique(.normalise_adduct_name(adducts))
        available <- .normalise_adduct_name(rules$name)
        missing_rules <- setdiff(requested, available)
        if (length(missing_rules)) {
            stop(
                "Unknown or unavailable adduct(s) for ", polarity, " mode: ",
                paste(missing_rules, collapse = ", ")
            )
        }
        rules <- rules[match(requested, available), , drop = FALSE]
    }
    if (!nrow(rules)) {
        stop("No adduct rules remain after filtering.")
    }

    compounds <- .normalise_metabolite_db(
        db,
        collapseIsomers = collapseIsomers,
        inferStructure = inferStructure, structureBackend = structureBackend,
        BPPARAM = BPPARAM,
        retainStructureRecords = any(rules$water_loss > 0L)
    )
    structure_records <- attr(compounds, "structure_records", exact = TRUE)
    structure_record_rows <- attr(compounds, "structure_record_rows", exact = TRUE)
    attr(compounds, "structure_records") <- NULL
    attr(compounds, "structure_record_rows") <- NULL
    if (!nrow(compounds)) {
        stop("No valid formula/monoisotopic-mass records were found in db.")
    }

    reactive_columns <- unique(rules$site_count_column[
        rules$min_reactive_sites > 0L & !is.na(rules$site_count_column)
    ])
    missing_reactive_columns <- reactive_columns[vapply(
        reactive_columns,
        function(column) {
            !column %in% names(compounds) || !any(is.finite(compounds[[column]]))
        },
        logical(1)
    )]
    if (length(missing_reactive_columns)) {
        warning(
            "No usable structure-derived reaction-site counts were found in: ",
            paste(missing_reactive_columns, collapse = ", "), ". Reactive-product ",
            "candidates will be retained with status 'unknown' and downweighted; ",
            "add these columns to db for chemical eligibility pruning.",
            call. = FALSE
        )
    }

    expected_parts <- vector("list", nrow(rules))
    compound_parts <- vector("list", nrow(rules))
    rule_parts <- vector("list", nrow(rules))
    max_h <- compounds$max_exchangeable_protons

    for (j in seq_len(nrow(rules))) {
        rule <- rules[j, , drop = FALSE]
        expected <- (rule$n_molecules * compounds$exactmass + rule$mass_shift) /
            abs(rule$charge)
        valid <- is.finite(expected) & expected > 0
        have_bound <- !is.na(max_h)
        hard_bound <- have_bound & (
            compounds$proton_bound_source == "provided" | max_h > 0
        )
        valid[hard_bound] <- valid[hard_bound] &
            rule$loss_h <= rule$n_molecules * max_h[hard_bound]
        structure_prior <- .adduct_structure_prior(compounds, rule)
        structure_known <- structure_prior$status != "unavailable"
        valid[structure_known] <- valid[structure_known] &
            structure_prior$score[structure_known] >= minStructureScore
        required_sites <- rule$min_reactive_sites[[1]]
        site_column <- rule$site_count_column[[1]]
        if (required_sites > 0L && !is.na(site_column) &&
            site_column %in% names(compounds)) {
            site_count <- compounds[[site_column]]
            known_site_count <- is.finite(site_count)
            valid[known_site_count] <- valid[known_site_count] &
                site_count[known_site_count] >= required_sites
        }
        idx <- which(valid)
        expected_parts[[j]] <- expected[idx]
        compound_parts[[j]] <- idx
        rule_parts[[j]] <- rep.int(j, length(idx))
    }

    expected_mz <- unlist(expected_parts, use.names = FALSE)
    compound_idx <- as.integer(unlist(compound_parts, use.names = FALSE))
    rule_idx <- as.integer(unlist(rule_parts, use.names = FALSE))
    ord <- order(expected_mz, method = "radix")

    structure(
        list(
            expected_mz = expected_mz[ord],
            compound_idx = compound_idx[ord],
            rule_idx = rule_idx[ord],
            compounds = compounds,
            structure_records = structure_records,
            structure_record_rows = structure_record_rows,
            rules = rules,
            family_rules = family_rules,
            polarity = polarity,
            maldi_matrix = if (is.null(maldiMatrix)) {
                "unspecified"
            } else {
                .normalise_maldi_matrix(maldiMatrix)
            },
            infer_structure = inferStructure,
            structure_backend = structureBackend,
            min_structure_score = minStructureScore,
            database_metadata = attr(db, "spamtp_database", exact = TRUE)
        ),
        class = "spamtp_mz_index"
    )
}

#' Print an indexed metabolite annotation search space
#'
#' @param x A `spamtp_mz_index` returned by [buildMZAnnotationIndex()].
#' @param ... Additional arguments reserved for print-method compatibility.
#'
#' @return `x`, invisibly.
#' @method print spamtp_mz_index
#' @export
print.spamtp_mz_index <- function(x, ...) {
    cat("SpaMTP m/z annotation index\n")
    cat("  polarity: ", x$polarity, "\n", sep = "")
    matrix_profile <- if (is.null(x$maldi_matrix)) "unspecified" else x$maldi_matrix
    cat("  MALDI matrix profile: ", matrix_profile, "\n", sep = "")
    cat("  metabolites: ", nrow(x$compounds), "\n", sep = "")
    parsed <- sum(!is.na(x$compounds$structure_valid) &
        x$compounds$structure_valid)
    cat("  SMILES structures parsed: ", parsed, "\n", sep = "")
    cat("  ion/reaction rules: ", nrow(x$rules), "\n", sep = "")
    cat("  expected ions: ", length(x$expected_mz), "\n", sep = "")
    invisible(x)
}

.normalise_ms1_spectrum <- function(ms1Spectrum) {
    if (is.null(ms1Spectrum)) {
        return(NULL)
    }
    if (is.numeric(ms1Spectrum) && is.null(dim(ms1Spectrum))) {
        ms1Spectrum <- data.frame(mz = ms1Spectrum, intensity = rep(1, length(ms1Spectrum)))
    }
    spectrum <- as.data.frame(ms1Spectrum)
    mz_col <- .first_existing_column(
        spectrum, c("mz", "m.z", "mass", "raw_mz", "observed_mz")
    )
    intensity_col <- .first_existing_column(
        spectrum, c("intensity", "abundance", "counts", "signal"),
        required = FALSE
    )
    intensity <- if (is.null(intensity_col)) {
        rep(1, nrow(spectrum))
    } else {
        suppressWarnings(as.numeric(spectrum[[intensity_col]]))
    }
    out <- data.frame(
        mz = suppressWarnings(as.numeric(spectrum[[mz_col]])),
        intensity = intensity
    )
    out <- out[is.finite(out$mz) & out$mz > 0, , drop = FALSE]
    out[order(out$mz), , drop = FALSE]
}

.find_spectrum_peak <- function(spectrum, targetMz, ppm) {
    if (is.null(spectrum) || !nrow(spectrum) || !is.finite(targetMz)) {
        return(c(mz = NA_real_, intensity = NA_real_, ppm = NA_real_))
    }
    pos <- findInterval(targetMz, spectrum$mz)
    idx <- unique(pmax(1L, pmin(nrow(spectrum), c(pos, pos + 1L))))
    errors <- abs(spectrum$mz[idx] - targetMz) / targetMz * 1e6
    best <- idx[which.min(errors)]
    if (!length(best) || errors[which.min(errors)] > ppm) {
        return(c(mz = NA_real_, intensity = NA_real_, ppm = NA_real_))
    }
    c(
        mz = spectrum$mz[best], intensity = spectrum$intensity[best],
        ppm = errors[which.min(errors)]
    )
}

.formula_mass_score <- function(observedMz, formula, rule, cache = NULL, ppm = 5) {
    counts <- .formula_counts(formula, cache)
    if (is.null(counts)) return(NA_real_)
    if (!all(names(counts) %in% c("C", "H", "O"))) return(1)
    # Negative fractional defects are legitimate for oxygen-rich acids and
    # large lipids. Check the actual CHO composition, including the ion rule,
    # instead of treating a wrapped nominal mass as chemical evidence.
    m <- .spamtp_mass_constants()
    masses <- c(C = m$carbon, H = m$hydrogen, O = m$oxygen)
    neutral_mass <- sum(counts * masses[names(counts)])
    predicted <- (rule$n_molecules * neutral_mass + rule$mass_shift) / abs(rule$charge)
    error <- abs(observedMz - predicted) / predicted * 1e6
    # Do not impose sub-5-ppm formula accuracy on rounded database masses.
    if (is.finite(error) && error > max(5, ppm)) 0.2 else 1
}

.ratio_score <- function(observed, expected) {
    if (!is.finite(observed) || observed <= 0 || !is.finite(expected) || expected <= 0) {
        return(0.5)
    }
    max(0.05, exp(-abs(log(observed / expected))))
}

.isotope_score <- function(
    expectedMz, formula, rule, spectrum, ppm, cache = NULL
) {
    if (is.null(spectrum) || is.null(.formula_counts(formula, cache))) {
        return(NA_real_)
    }
    m <- .spamtp_mass_constants()
    charge <- abs(rule$charge)
    main <- .find_spectrum_peak(spectrum, expectedMz, ppm)
    score <- 1

    carbon_count <- rule$n_molecules * .element_count(formula, "C", cache) +
        (rule$carbon_delta %||% 0)
    if (carbon_count > 0) {
        isotope <- .find_spectrum_peak(
            spectrum, expectedMz + m$carbon13_spacing / charge, ppm
        )
        if (!is.finite(isotope[["mz"]])) {
            score <- score * 0.5
        } else if (is.finite(main[["intensity"]]) && main[["intensity"]] > 0) {
            ratio <- isotope[["intensity"]] / main[["intensity"]]
            score <- score * .ratio_score(ratio, carbon_count * 0.0107)
        }
    }

    n_cl <- rule$n_molecules * .element_count(formula, "Cl", cache) +
        (rule$chlorine_delta %||% as.integer(isTRUE(rule$name == "M+Cl")))
    n_br <- rule$n_molecules * .element_count(formula, "Br", cache) +
        (rule$bromine_delta %||% as.integer(isTRUE(rule$name == "M+Br")))
    if (n_cl > 0 || n_br > 0) {
        spacing <- if (n_br > 0) m$bromine_m2_spacing else m$chlorine_m2_spacing
        expected_ratio <- if (n_br > 0) n_br * 0.97 else n_cl * 0.32
        isotope_m2 <- .find_spectrum_peak(
            spectrum, expectedMz + spacing / charge, ppm
        )
        if (!is.finite(isotope_m2[["mz"]])) {
            score <- score * 0.2
        } else if (is.finite(main[["intensity"]]) && main[["intensity"]] > 0) {
            score <- score * .ratio_score(
                isotope_m2[["intensity"]] / main[["intensity"]], expected_ratio
            )
        }
    }
    score
}

.adduct_family_score <- function(neutralMass, rule, allRules, spectrum, ppm) {
    complex_types <- c("multimer", "solvent", "metal_exchange", "matrix_adduct")
    if (is.null(spectrum) || !rule$complexity %in% complex_types) {
        return(NA_real_)
    }
    base_names <- unique(c(
        rule$base_adduct,
        if (rule$polarity == "positive") c("M+H", "M+Na") else "M-H"
    ))
    base_rules <- allRules[allRules$name %in% base_names &
        allRules$n_molecules == 1, , drop = FALSE]
    if (!nrow(base_rules)) {
        return(NA_real_)
    }
    base_mz <- (base_rules$n_molecules * neutralMass + base_rules$mass_shift) /
        abs(base_rules$charge)
    found <- vapply(base_mz, function(x) {
        is.finite(.find_spectrum_peak(spectrum, x, ppm)[["mz"]])
    }, logical(1))
    if (any(found)) 1 else 0.1
}

.empty_annotation_result <- function() {
    data.frame(
        observed_mz = numeric(), expected_mz = numeric(), ppm_error = numeric(),
        score = numeric(), mass_score = numeric(), rule_prior = numeric(),
        water_loss = integer(), dehydration_sites = numeric(),
        max_water_losses = numeric(), dehydration_status = character(),
        chemical_score = numeric(), structure_score = numeric(),
        structure_status = character(), structure_rule_evidence = character(),
        structure_evidence = character(),
        positive_site_atoms = character(), negative_site_atoms = character(),
        alkali_site_atoms = character(),
        reactive_site_score = numeric(), reactive_site_status = character(),
        isotope_score = numeric(),
        adduct_network_score = numeric(), adduct = character(), charge = integer(),
        maldi_matrix = character(), rule_class = character(),
        rule_source = character(), reactive_group = character(),
        min_reactive_sites = integer(),
        matrix_target_groups = character(), matrix_target_sites = numeric(),
        matrix_target_status = character(),
        formula = character(), formula_status = character(),
        neutral_mass = numeric(), metabolite_ids = character(),
        inchikeys = character(), metabolite_names = character(), ramp_ids = character(),
        stringsAsFactors = FALSE
    )
}

#' Query an indexed metabolite annotation search space
#'
#' Candidates are pruned by ppm tolerance and proton/charge validity during
#' index construction, then ranked using mass error, rule prior, formula/ion mass consistency,
#' isotope, and contextual adduct-family evidence. Positive water-loss
#' candidates additionally require compatible alcohol OH and neighbouring
#' carbon H sites. Older structure tables are checked lazily for mass-matched
#' records; `inferStructure = "never"` requires precomputed dehydration
#' capacity. Isomer identities are resolved separately for these losses.
#'
#' @param observedMz Numeric vector of observed m/z values.
#' @param index A `spamtp_mz_index` from [buildMZAnnotationIndex()].
#' @param ppm Mass tolerance in parts per million.
#' @param ms1Spectrum Optional contextual spectrum with `mz` and `intensity`
#'   columns. It should represent the same retention-time window or spatial
#'   pixel/region as the queried peaks.
#' @param checkFormulaMass Check CHO formula/ion mass consistency. A mismatch
#'   beyond the query tolerance (at least 5 ppm for rounded reference masses)
#'   receives a 0.2 chemical score. Negative fractional mass defects alone
#'   are not penalized.
#' @param checkIsotopes Score carbon-13 and, when relevant, chlorine/bromine
#'   M+2 patterns.
#' @param checkAdductNetwork Downweight complex ions whose base monomer
#'   family is absent from `ms1Spectrum`.
#' @param minScore Minimum final score to retain.
#'
#' @return A ranked candidate data frame. `water_loss` is the required number
#'   of neutral water losses; `dehydration_sites` and `max_water_losses` record
#'   the checked structure's site count and capacity. `dehydration_status` is
#'   `"structure-compatible"` for retained positive water losses and
#'   `"not_required"` for other rules. Missing capacities for other rules are
#'   left as `NA`. These fields describe eligibility, not confirmed identity.
#'   `formula_status` is `"parsed"`, `"unavailable"`, or `"not_checked"`.
#'   Unsupported formulas retain exact-mass candidates, with unavailable
#'   formula-based scores recorded as `NA` and omitted from the score product.
#'   Water-loss candidates require complete formula counts as well as structure
#'   evidence. A single `spamtp_formula_warning` summarizes unavailable formulas
#'   per query; its `formulas` field contains the affected reference strings.
#' @examples
#' data(spamtp_example)
#' index <- buildMZAnnotationIndex(spamtp_example$metabolites, adducts = "M+H")
#' queryMZAnnotationIndex(c(149.04445, 181.07066), index, ppm = 5)
#' @export
queryMZAnnotationIndex <- function(
    observedMz, index, ppm = 5,
    ms1Spectrum = NULL,
    checkFormulaMass = TRUE,
    checkIsotopes = TRUE,
    checkAdductNetwork = TRUE,
    minScore = 0
) {
    if (!inherits(index, "spamtp_mz_index")) {
        stop("index must be created by buildMZAnnotationIndex().")
    }
    observedMz <- suppressWarnings(as.numeric(observedMz))
    if (!length(observedMz) || any(!is.finite(observedMz) | observedMz <= 0)) {
        stop("observedMz must contain positive finite values.")
    }
    if (length(ppm) != 1 || !is.finite(ppm) || ppm < 0) {
        stop("ppm must be a single non-negative number.")
    }
    spectrum <- .normalise_ms1_spectrum(ms1Spectrum)
    tolerance <- ppm / 1e6
    results <- vector("list", length(observedMz))
    dehydration_cache <- new.env(parent = emptyenv())
    formula_cache <- new.env(parent = emptyenv())
    lower <- observedMz / (1 + tolerance)
    upper <- if (tolerance < 1) observedMz / (1 - tolerance) else rep(Inf, length(lower))
    lower_open <- lower - .Machine$double.eps * pmax(1, abs(lower))
    lower_index <- findInterval(lower_open, index$expected_mz)
    upper_index <- findInterval(upper, index$expected_mz)

    for (i in seq_along(observedMz)) {
        obs <- observedMz[[i]]
        lo <- lower_index[[i]]
        hi <- upper_index[[i]]
        if (hi <= lo) next
        hit <- seq.int(lo + 1L, hi)
        expected <- index$expected_mz[hit]
        error <- abs(obs - expected) / expected * 1e6
        keep <- is.finite(error) & error <= ppm
        if (!any(keep)) next
        hit <- hit[keep]
        expected <- expected[keep]
        error <- error[keep]
        compounds <- index$compounds[index$compound_idx[hit], , drop = FALSE]
        rules <- index$rules[index$rule_idx[hit], , drop = FALSE]
        if (is.null(rules$water_loss)) rules$water_loss <- rep(0L, nrow(rules))
        if (any(rules$water_loss > 0L) && !is.null(index$structure_records)) {
            members <- lapply(seq_along(hit), function(j) {
                if (rules$water_loss[[j]] > 0L) {
                    index$structure_record_rows[[index$compound_idx[hit[[j]]]]]
                } else {
                    NA_integer_
                }
            })
            source_rows <- rep(seq_along(hit), lengths(members))
            record_rows <- unlist(members, use.names = FALSE)
            compounds <- compounds[source_rows, , drop = FALSE]
            original <- !is.na(record_rows)
            compounds[original, ] <- index$structure_records[
                record_rows[original], names(compounds),
                drop = FALSE
            ]
            rules <- rules[source_rows, , drop = FALSE]
            expected <- expected[source_rows]
            error <- error[source_rows]
        }
        dehydration <- .resolve_dehydration(
            compounds, rules, index$infer_structure %||% "auto",
            dehydration_cache, formula_cache
        )
        keep <- dehydration$keep
        if (!any(keep)) next
        compounds <- dehydration$compounds[keep, , drop = FALSE]
        dehydration_status <- dehydration$status[keep]
        rules <- rules[keep, , drop = FALSE]
        expected <- expected[keep]
        error <- error[keep]
        structure_prior <- lapply(seq_along(error), function(j) {
            .adduct_structure_prior(
                compounds[j, , drop = FALSE], rules[j, , drop = FALSE]
            )
        })
        structure_score <- vapply(
            structure_prior, function(x) x$score[[1]], numeric(1)
        )
        structure_status <- vapply(
            structure_prior, function(x) x$status[[1]], character(1)
        )
        structure_rule_evidence <- vapply(
            structure_prior, function(x) x$rationale[[1]], character(1)
        )
        matrix_compatibility <- .matrix_structure_compatibility(
            compounds, index$maldi_matrix
        )

        mass_score <- if (ppm == 0) {
            rep(1, length(error))
        } else {
            exp(-0.5 * (error / (ppm / 3))^2)
        }
        chemical_score <- if (isTRUE(checkFormulaMass)) {
            vapply(seq_along(error), function(j) {
                .formula_mass_score(
                    obs,
                    compounds$formula[[j]],
                    rules[j, , drop = FALSE], formula_cache, ppm = ppm
                )
            }, numeric(1))
        } else {
            rep(1, length(error))
        }
        site_status <- vapply(seq_along(error), function(j) {
            required_sites <- rules$min_reactive_sites[[j]]
            site_column <- rules$site_count_column[[j]]
            if (required_sites <= 0L) {
                return("not_required")
            }
            if (is.na(site_column) || !site_column %in% names(compounds)) {
                return("unknown")
            }
            count <- compounds[[site_column]][[j]]
            if (!is.finite(count)) {
                return("unknown")
            }
            if (count >= required_sites) "verified" else "insufficient"
        }, character(1))
        reactive_site_score <- ifelse(site_status == "unknown", 0.25, 1)
        isotope_score <- if (isTRUE(checkIsotopes)) {
            vapply(seq_along(error), function(j) {
                .isotope_score(
                    expected[[j]], compounds$formula[[j]],
                    rules[j, , drop = FALSE], spectrum, ppm, formula_cache
                )
            }, numeric(1))
        } else {
            rep(NA_real_, length(error))
        }
        network_score <- if (isTRUE(checkAdductNetwork)) {
            vapply(seq_along(error), function(j) {
                .adduct_family_score(
                    compounds$exactmass[[j]], rules[j, , drop = FALSE],
                    index$family_rules, spectrum, ppm
                )
            }, numeric(1))
        } else {
            rep(NA_real_, length(error))
        }

        formula_checked <- rules$water_loss > 0L |
            isTRUE(checkFormulaMass) |
            (isTRUE(checkIsotopes) && !is.null(spectrum))
        formula_status <- vapply(seq_along(error), function(j) {
            if (!formula_checked[[j]]) {
                return("not_checked")
            }
            if (is.null(.formula_counts(compounds$formula[[j]], formula_cache))) {
                "unavailable"
            } else {
                "parsed"
            }
        }, character(1))
        evidence_chemical <- ifelse(is.na(chemical_score), 1, chemical_score)
        evidence_isotope <- ifelse(is.na(isotope_score), 1, isotope_score)
        evidence_network <- ifelse(is.na(network_score), 1, network_score)
        final_score <- mass_score * rules$prior * evidence_chemical * structure_score *
            reactive_site_score *
            evidence_isotope * evidence_network

        out <- data.frame(
            observed_mz = obs,
            expected_mz = expected,
            ppm_error = error,
            score = final_score,
            mass_score = mass_score,
            rule_prior = rules$prior,
            water_loss = rules$water_loss,
            dehydration_sites = compounds$dehydration_sites,
            max_water_losses = compounds$max_water_losses,
            dehydration_status = dehydration_status,
            chemical_score = chemical_score,
            structure_score = structure_score,
            structure_status = structure_status,
            structure_rule_evidence = structure_rule_evidence,
            structure_evidence = compounds$structure_evidence,
            positive_site_atoms = compounds$positive_site_atoms,
            negative_site_atoms = compounds$negative_site_atoms,
            alkali_site_atoms = compounds$alkali_site_atoms,
            reactive_site_score = reactive_site_score,
            reactive_site_status = site_status,
            isotope_score = isotope_score,
            adduct_network_score = network_score,
            adduct = rules$name,
            charge = rules$charge,
            maldi_matrix = rules$maldi_matrix,
            rule_class = rules$rule_class,
            rule_source = rules$rule_source,
            reactive_group = rules$reactive_group,
            min_reactive_sites = rules$min_reactive_sites,
            matrix_target_groups = matrix_compatibility$target,
            matrix_target_sites = matrix_compatibility$sites,
            matrix_target_status = matrix_compatibility$status,
            formula = compounds$formula,
            formula_status = unname(formula_status),
            neutral_mass = compounds$exactmass,
            metabolite_ids = compounds$isomers,
            inchikeys = compounds$isomers_inchikey,
            metabolite_names = compounds$isomers_names,
            ramp_ids = compounds$ramp_ids,
            stringsAsFactors = FALSE
        )
        allowed <- out$score >= minScore
        # Apply the threshold again after lazy structure resolution.
        water <- out$water_loss > 0L
        allowed[water] <- allowed[water] &
            out$structure_score[water] >= (index$min_structure_score %||% 0.05)
        out <- out[allowed, , drop = FALSE]
        duplicate_loss <- out$water_loss > 0L & duplicated(out)
        results[[i]] <- out[!duplicate_loss, , drop = FALSE]
    }

    unavailable <- vapply(as.list(formula_cache), is.null, logical(1))
    .warn_formula_evidence(sub("^formula:", "", names(unavailable)[unavailable]))
    result <- do.call(rbind, results)
    if (is.null(result) || !nrow(result)) {
        return(.empty_annotation_result())
    }
    result <- result[order(result$observed_mz, -result$score, result$ppm_error), , drop = FALSE]
    rownames(result) <- NULL
    result
}

# Reusing an index can narrow its rule space, but cannot rebuild its database
# or structural evidence. Reject conflicting construction options explicitly.
.reuse_mz_index <- function(
    index, adducts = NULL, db = NULL, rules = NULL,
    inferStructure = NULL, structureBackend = NULL,
    minStructureScore = NULL, suppliedBpparam = FALSE
) {
    if (!is.null(db) || !is.null(rules) || suppliedBpparam) {
        stop("With index supplied, omit db, rules and BPPARAM; use buildMZAnnotationIndex() to change them.",
            call. = FALSE)
    }
    requested <- list(infer_structure = inferStructure,
        structure_backend = structureBackend, min_structure_score = minStructureScore)
    argument_names <- c(infer_structure = "inferStructure",
        structure_backend = "structureBackend", min_structure_score = "minStructureScore")
    for (name in names(requested)) {
        value <- requested[[name]]
        if (!is.null(value) && !identical(value, index[[name]])) {
            stop("The supplied index does not match ", argument_names[[name]],
                "; rebuild it with buildMZAnnotationIndex().", call. = FALSE)
        }
    }
    if (!is.null(adducts)) {
        requested <- unique(.normalise_adduct_name(adducts))
        available <- .normalise_adduct_name(index$rules$name)
        missing_rules <- setdiff(requested, available)
        if (!length(requested) || length(missing_rules)) {
            stop("Requested adduct(s) are unavailable in the supplied index: ",
                paste(missing_rules, collapse = ", "), ". Rebuild the index to add rules.",
                call. = FALSE)
        }
        selected <- match(requested, available)
        keep <- index$rule_idx %in% selected
        index$expected_mz <- index$expected_mz[keep]
        index$compound_idx <- index$compound_idx[keep]
        index$rule_idx <- match(index$rule_idx[keep], selected)
        index$rules <- index$rules[selected, , drop = FALSE]
    }
    index
}

#' Annotate one or more observed m/z values
#'
#' Convenience wrapper that builds an index when one is not supplied. Reuse a
#' pre-built index for repeated or large searches.
#'
#' @param observedMz Numeric vector of observed m/z values.
#' @param db Metabolite database used when `index` is `NULL`.
#' @param index Optional pre-built `spamtp_mz_index`. `adducts` can restrict
#'   the existing rules; unavailable rules or conflicting construction options
#'   require a new index. Omit `db`, `rules` and `BPPARAM` when reusing an index.
#' @param polarity Ion mode. When `NULL`, infer it from `index` or the matrix
#'   profile, falling back to positive mode.
#' @param adducts Optional adduct subset. When `NULL`, retain the complete rule
#'   space selected by the matrix profile or polarity.
#' @param rules Optional custom rule table.
#' @param maldiMatrix Optional MALDI matrix/reagent profile used to select
#'   rules automatically when `rules` is `NULL`.
#' @param ppm Mass tolerance in ppm.
#' @param ms1Spectrum Optional contextual spectrum.
#' @param databaseVersion SpaMTPdb/RaMP version used when `db = NULL`.
#' @param databaseSource Database source used when `db = NULL`; see
#'   [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @param inferStructure,structureBackend,BPPARAM,min_structure_score
#'   Structure-aware
#'   rule-selection arguments passed to [buildMZAnnotationIndex()].
#' @param ... Additional arguments passed to [queryMZAnnotationIndex()].
#'
#' @return A ranked candidate data frame.
#' @examples
#' data(spamtp_example)
#' index <- buildMZAnnotationIndex(spamtp_example$metabolites, adducts = "M+H")
#' annotateMZ(149.04445, index = index, ppm = 5)
#' @export
annotateMZ <- function(
    observedMz, db = NULL, index = NULL,
    polarity = NULL,
    adducts = NULL, rules = NULL, maldiMatrix = NULL,
    ppm = 5,
    ms1Spectrum = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    inferStructure = c("auto", "never", "always"),
    structureBackend = c("auto", "native"),
    BPPARAM = BiocParallel::SerialParam(),
    minStructureScore = 0.05, ...
) {
    if (!is.null(index) && !inherits(index, "spamtp_mz_index")) {
        stop("index must be created by buildMZAnnotationIndex().")
    }
    polarity <- if (is.null(polarity) && inherits(index, "spamtp_mz_index")) {
        index$polarity
    } else {
        .resolve_maldi_polarity(polarity, maldiMatrix)
    }
    if (is.null(index)) {
        if (is.null(db)) {
            db <- .spamtp_db_resource(
                "chem_props",
                version = databaseVersion,
                source = match.arg(databaseSource),
                localDir = databaseLocalDir
            )
        }
        index <- buildMZAnnotationIndex(
            db = db, polarity = polarity, adducts = adducts, rules = rules,
            maldiMatrix = maldiMatrix, inferStructure = inferStructure,
            structureBackend = structureBackend,
            BPPARAM = BPPARAM,
            minStructureScore = minStructureScore
        )
    } else if (!identical(index$polarity, polarity)) {
        stop("The supplied index polarity does not match polarity.")
    } else if (!is.null(maldiMatrix) &&
        !identical(index$maldi_matrix, .normalise_maldi_matrix(maldiMatrix))) {
        stop("The supplied index MALDI matrix profile does not match maldi_matrix.")
    } else {
        index <- .reuse_mz_index(index, adducts, db, rules,
            inferStructure = if (!missing(inferStructure)) match.arg(inferStructure),
            structureBackend = if (!missing(structureBackend)) match.arg(structureBackend),
            minStructureScore = if (!missing(minStructureScore)) minStructureScore,
            suppliedBpparam = !missing(BPPARAM))
    }
    queryMZAnnotationIndex(
        observedMz = observedMz, index = index, ppm = ppm,
        ms1Spectrum = ms1Spectrum, ...
    )
}
