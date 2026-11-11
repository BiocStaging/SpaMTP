#### SpaMTP m/z Annotation Functions #####################################################################################################################################################################################

#' Subset a Bioconductor experiment by a list of m/z's
#'
#' @param data A Bioconductor experiment for subsetting.
#' @param features A list of character strings defining the features/mz values
#'   to subset against.
#' @param assay A character string identifying the primary or alternative
#'   experiment which contains the count data being subset.
#'
#' @returns A subset Bioconductor experiment containing only m/z values that
#'   were specified
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' subsetMZFeatures(spe, c("ion1", "ion3"))
subsetMZFeatures <- function(data, features, assay = "main") {
    experiment <- .experimentForAssay(data, assay)
    features <- as.character(features)
    if (anyNA(features) || anyDuplicated(features) ||
        !all(features %in% rownames(experiment))) {
        stop(
            "features must be unique names in the selected experiment.",
            call. = FALSE
        )
    }
    .replaceExperiment(data, experiment[features, , drop = FALSE], assay)
}


#' Filters the annotation list to only include the first n number of annotations
#' per m/z
#'
#' @param annotationColumn Vector from the `rowData()` column containing m/z
#'   annotations.
#' @param n Numeric value defining the number of annotations to keep (default =
#'   3).
#'
#' @return Vector containing the first n number of annotations
#'
#' @noRd
labels_to_show <- function(annotationColumn, n = 3) {
    new_column <- vapply(strsplit(annotationColumn, "; "), function(x) {
        # Filter out NA values and select the first three entries
        non_na_values <- x[!is.na(x)]
        paste(utils::head(non_na_values, n), collapse = "; ")
    }, character(1))
    return(new_column)
}


#' Annotates m/z values stored in a SpaMTP Object
#'
#' This function assigns each valid m/z peak with one/multiple metabolite names
#' based on the mass difference between the observed value and the theoretical
#' value documented in the reference database.
#' Versioned RaMP, HMDB, LIPID MAPS, ChEBI, and GNPS resources are provided by
#' SpaMTPdb and can be combined with user-supplied reference tables.
#'
#' @param data Bioconductor experiment containing m/z values for annotation.
#' @param db Reference metabolite dataset in the form of a data.frame. When
#'   `NULL`, the versioned SpaMTPdb `chem_props` table is used, unless a
#'   pre-built `index` is supplied through `...`. Versioned resources are loaded
#'   from SpaMTPdb and may be staged locally for offline use.
#' @param assay Primary or alternative experiment containing the m/z counts
#'   being annotated. Defaults to `"main"`.
#' @param rawMzColumn Character string naming the feature-metadata column
#'   containing raw m/z values without the `mz-` prefix. This is created by
#'   `cardinalToSpatialExperiment()` (default = "raw_mz"). If omitted, mass
#'   metadata or numeric feature names are used as a fallback. An explicitly
#'   named column must exist. Feature rownames must be unique and non-empty.
#' @param ppm Mass tolerance in ppm (default 5). Set to zero
#'   for exact numerical matches.
#' @param adducts Optional adduct names/notations; see `adductRules()`. If
#'   `NULL`, use the complete rule space selected from `maldiMatrix`, or all
#'   validated general rules for the selected polarity when no matrix is given.
#' @param polarity Character string defining the ion mode. When `NULL`, use the
#'   MALDI matrix profile default, a supplied index's mode, or positive mode
#'   when neither is available.
#' @param filepath Optional output filename prefix, including its directory.
#'   For example, `file.path("results", "sample")` writes
#'   `results/sample_annotated_mz_peaks.csv`. NULL disables export.
#' @param returnOnlyAnnotated Boolean value indicating if the annotated
#'   Bioconductor experiment should only include m/z values that were
#'   successfully annotated (default = TRUE).
#' @param saveIntermediate Boolean indicating whether to store the scored
#'   annotation result and its pipeline/RaMP provenance in
#'   `S4Vectors::metadata(data)$mz_annotation`. A compatibility copy
#'   is retained under the `db_3` key (default = TRUE).
#' @param minScore Minimum annotation score retained in the stored candidate
#'   table. The default `0` keeps all ppm-valid candidates so downstream
#'   pathway functions can apply a user-defined threshold without re-running
#'   annotation.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param maldiMatrix Optional MALDI matrix or derivatization reagent name.
#'   When supplied, SpaMTP selects validated matrix-specific rules
#'   automatically. `adducts` is optional and only restricts that automatic
#'   search space when explicitly supplied.
#' @param databaseVersion SpaMTPdb/RaMP version used when `db = NULL`.
#' @param databaseSource Database source used when `db = NULL`; see
#'   [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @param ... Additional indexed annotation/scoring arguments passed to
#'   `annotateTable()`, such as `index`, `rules`, `ms1_spectrum`,
#'   `infer_structure`, `structure_backend`, `BPPARAM`, or
#'   `min_structure_score`.
#'
#' @returns A Bioconductor experiment with m/z values annotated. These
#'   annotations are stored in the selected experiment's rowData().
#'   Stored candidates include `FormulaStatus`: `"parsed"`, `"unavailable"`,
#'   or `"not_checked"`. One warning summarizes unavailable formula-based
#'   evidence; inspect these records with [annotationCandidates()].
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' index <- buildMZAnnotationIndex(spamtp_example$metabolites, adducts = "M+H")
#' annotated <- annotateSM(spe,
#'     index = index, rawMzColumn = "mz", ppm = 5,
#'     returnOnlyAnnotated = FALSE, verbose = FALSE
#' )
#' head(annotationCandidates(annotated))
annotateSM <- function(
    data,
    db = NULL,
    assay = "main",
    rawMzColumn = "raw_mz",
    ppm = 5,
    adducts = NULL,
    polarity = NULL,
    filepath = NULL,
    returnOnlyAnnotated = TRUE,
    saveIntermediate = TRUE,
    minScore = 0,
    verbose = TRUE,
    maldiMatrix = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    ...
) {
    experiment <- .experimentForAssay(data, assay)
    if (is.null(rownames(experiment)) || anyNA(rownames(experiment)) ||
        any(!nzchar(rownames(experiment))) || anyDuplicated(rownames(experiment))) {
        stop("Annotation requires unique, non-empty feature rownames; set rownames() before annotateSM().",
            call. = FALSE)
    }

    annotation_args <- list(...)
    supplied_index <- annotation_args$index
    polarity <- if (is.null(polarity) &&
        inherits(supplied_index, "spamtp_mz_index")) {
        supplied_index$polarity
    } else {
        .resolve_maldi_polarity(polarity, maldiMatrix)
    }
    resolved_matrix <- if (!is.null(maldiMatrix)) {
        .normalise_maldi_matrix(maldiMatrix)
    } else if (inherits(supplied_index, "spamtp_mz_index")) {
        supplied_index$maldi_matrix
    } else {
        "unspecified"
    }
    default_ramp <- is.null(db) && is.null(supplied_index)
    if (default_ramp) {
        db <- .spamtp_db_resource(
            "chem_props",
            version = databaseVersion,
            source = match.arg(databaseSource),
            localDir = databaseLocalDir
        )
        verbose_message(
            messageText = paste0(
                "Using ",
                .spamtp_db_label(db),
                " for annotation ... "
            ),
            verbose = verbose
        )
    }

    ## Extract m/z values from assay feature metadata.
    featureMetadata <- .featureMetadata(data, assay)
    preserve_ids <- saveIntermediate ||
        ".spamtp_feature_id" %in% names(featureMetadata)
    featureMetadata$.spamtp_feature_id <- .annotationFeatureIds(
        featureMetadata,
        initialize = TRUE
    )
    mzValues <- featureMetadata[[rawMzColumn]]
    if (is.null(mzValues)) {
        if (!missing(rawMzColumn)) {
            stop("Feature metadata column '", rawMzColumn, "' is absent from assay '", assay, "'.",
                call. = FALSE)
        }
        mzValues <- .massValues(featureMetadata, rownames(experiment))
    }
    if (!is.numeric(mzValues) || length(mzValues) != nrow(experiment) ||
        any(!is.finite(mzValues))) {
        stop(
            "Annotation requires one finite numeric m/z value per feature.",
            call. = FALSE
        )
    }
    mz_df <- data.frame(row_id = seq_along(mzValues), mz = mzValues)

    db_3 <- annotateTable(
        mzDf = mz_df,
        db = db,
        ppm = ppm,
        adducts = adducts,
        polarity = polarity,
        verbose = verbose,
        minScore = minScore,
        maldiMatrix = maldiMatrix,
        ...
    )
    # Keep identity separate from the display name, which users may rename.
    db_3$feature <- rownames(experiment)[db_3$ID]
    db_3$feature_id <- featureMetadata$.spamtp_feature_id[db_3$ID]

    if (saveIntermediate) {
        alternatives <- if (inherits(data, "SingleCellExperiment")) SingleCellExperiment::altExpNames(data) else character()
        stores <- .storedData(data, "modality_annotations") %||% list()
        previous <- .storedData(data, "mz_annotation")
        previous_legacy <- .storedData(data, "db_3")
        if (!is.null(previous) || !is.null(previous_legacy)) {
            previous_assay <- previous$metadata$assay %||% "main"
            previous_key <- if (previous_assay %in% alternatives) previous_assay else "main"
            stores[[previous_key]] <- list(mz_annotation = previous, db_3 = previous_legacy)
        }
        data <- .setStoredData(data, "mz_annotation", .annotation_store(
            db_3,
            metadata = list(
                assay = assay,
                raw_mz_column = rawMzColumn,
                polarity = polarity,
                maldi_matrix = resolved_matrix,
                ppm = ppm,
                min_score = minScore,
                ramp_version = .annotation_ramp_version(supplied_index %||% db),
                adducts = if (is.null(adducts)) {
                    if (inherits(supplied_index, "spamtp_mz_index")) {
                        "complete pre-built index rule space"
                    } else if (is.null(maldiMatrix)) {
                        "all validated general rules"
                    } else {
                        "automatically selected from MALDI matrix profile"
                    }
                } else {
                    adducts
                },
                database = if (default_ramp) {
                    .spamtp_db_label(db)
                } else if (!is.null(supplied_index)) {
                    "pre-built annotation index"
                } else {
                    "user-supplied database"
                }
            )
        ))
        # Keep the historical location for packages and serialized objects that
        # still inspect it directly. New pathway code reads mz_annotation first.
        data <- .setStoredData(data, "db_3", db_3)
        key <- if (assay %in% alternatives) assay else "main"
        stores[[key]] <- list(
            mz_annotation = .storedData(data, "mz_annotation"),
            db_3 = db_3
        )
        data <- .setStoredData(data, "modality_annotations", stores)
    }

    if (!(is.null(filepath))) {
        path <- paste0(filepath, "_annotated_mz_peaks.csv")
        message("Writing annotated m/z file to: '", path, "'")
        data.table::fwrite(db_3, path)
    }

    verbose_message(
        messageText = "Adding annotations to rowData .... ",
        verbose = verbose
    )

    result_df <- db_3 %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::summarise(
            all_IsomerNames = paste(.data$IsomerNames, collapse = "; "),
            all_Isomers = paste(.data$Isomers, collapse = "; "),
            all_Isomers_IDs = paste(.data$Isomers_IDs, collapse = "; "),
            all_Adducts = paste(unique(.data$Adduct), collapse = "; "),
            all_Formulas = paste(unique(.data$Formula), collapse = "; "),
            all_Errors = paste(round(.data$Error, 4), collapse = "; "),
            all_Scores = paste(round(.data$Score, 4), collapse = "; "),
            all_Ramp_IDs = paste(
                unique(.data$Ramp_IDs[!is.na(.data$Ramp_IDs)]),
                collapse = "; "
            )
        )
    result_df <- data.frame(result_df)
    matched <- match(featureMetadata$.spamtp_feature_id, result_df$feature_id)
    featureMetadata$mz_names <- rownames(experiment)
    for (column in setdiff(colnames(result_df), "feature_id")) {
        value <- result_df[[column]][matched]
        value[is.na(value)] <- "No Annotation"
        featureMetadata[[column]] <- value
    }
    if (!preserve_ids) featureMetadata$.spamtp_feature_id <- NULL
    data <- .setFeatureMetadata(data, featureMetadata, assay)
    if (isTRUE(returnOnlyAnnotated)) {
        data <- subsetMZFeatures(
            data,
            rownames(experiment)[!is.na(matched)],
            assay
        )
    }

    return(data)
}


#' @noRd
annotateTable <- function(
    mzDf, db = NULL, ppm = 5, adducts = NULL,
    polarity = NULL,
    verbose = TRUE, index = NULL, rules = NULL,
    ms1Spectrum = NULL, checkFormulaMass = TRUE,
    checkIsotopes = TRUE,
    checkAdductNetwork = TRUE, minScore = 0,
    maldiMatrix = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL,
    inferStructure = c("auto", "never", "always"),
    structureBackend = c("auto", "native"),
    BPPARAM = BiocParallel::SerialParam(),
    minStructureScore = 0.05
) {
    if (!is.data.frame(mzDf) || !all(c("row_id", "mz") %in% names(mzDf))) {
        stop("mzDf must be a data.frame containing row_id and mz columns.")
    }
    mzDf$mz <- suppressWarnings(as.numeric(mzDf$mz))
    if (anyNA(mzDf$row_id) || anyDuplicated(mzDf$row_id)) {
        stop("mzDf$row_id must contain unique, non-missing feature identifiers.", call. = FALSE)
    }
    if (any(!is.finite(mzDf$mz) | mzDf$mz <= 0)) {
        stop("mzDf$mz must contain positive finite values.")
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
        verbose_message(
            messageText = "Building sorted adduct candidate index ... ",
            verbose = verbose
        )
        index <- buildMZAnnotationIndex(
            db = db, polarity = polarity, adducts = adducts, rules = rules,
            maldiMatrix = maldiMatrix, inferStructure = inferStructure,
            structureBackend = structureBackend,
            BPPARAM = BPPARAM,
            minStructureScore = minStructureScore
        )
    } else if (!inherits(index, "spamtp_mz_index")) {
        stop("index must be created by buildMZAnnotationIndex().")
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

    verbose_message(
        messageText = paste0(
            "Searching candidate index at ",
            ppm,
            " ppm ... "
        ),
        verbose = verbose
    )
    candidates <- if (!nrow(mzDf)) .empty_annotation_result() else queryMZAnnotationIndex(
        observedMz = unique(mzDf$mz),
        index = index,
        ppm = ppm,
        ms1Spectrum = ms1Spectrum,
        checkFormulaMass = checkFormulaMass,
        checkIsotopes = checkIsotopes,
        checkAdductNetwork = checkAdductNetwork,
        minScore = minScore
    )

    empty_result <- function() {
        data.frame(
            ID = mzDf$row_id[0], Match = logical(), observed_mz = numeric(),
            Reference_mz = numeric(), Error = numeric(), Adduct = character(),
            Formula = character(), FormulaStatus = character(),
            Exactmass = numeric(), Isomers = character(),
            InchiKeys = character(), IsomerNames = character(),
            Isomers_IDs = character(), Ramp_IDs = character(), Score = numeric(),
            MassScore = numeric(), ChemicalScore = numeric(),
            WaterLoss = integer(), DehydrationSites = numeric(),
            MaxWaterLosses = numeric(), DehydrationStatus = character(),
            StructureScore = numeric(), StructureStatus = character(),
            StructureRuleEvidence = character(), StructureEvidence = character(),
            PositiveSiteAtoms = character(), NegativeSiteAtoms = character(),
            AlkaliSiteAtoms = character(),
            ReactiveSiteScore = numeric(), ReactiveSiteStatus = character(),
            IsotopeScore = numeric(), AdductNetworkScore = numeric(),
            MALDIMatrix = character(), RuleClass = character(),
            RuleSource = character(), ReactiveGroup = character(),
            MinimumReactiveSites = integer(),
            MatrixTargetGroups = character(), MatrixTargetSites = numeric(),
            MatrixTargetStatus = character(),
            stringsAsFactors = FALSE
        )
    }
    if (!nrow(candidates)) {
        return(empty_result())
    }

    row_map <- data.frame(ID = mzDf$row_id, observed_mz = mzDf$mz)
    candidates$.candidate_order <- seq_len(nrow(candidates))
    matched <- merge(
        row_map, candidates,
        by = "observed_mz", all = FALSE, sort = FALSE
    )
    matched <- matched[order(matched$ID, matched$.candidate_order), , drop = FALSE]

    standardise_ids <- function(value) {
        if (length(value) != 1 || is.na(value) || !nzchar(value)) {
            return("")
        }
        ids <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
        ids <- ids[nzchar(ids)]
        ids <- vapply(ids, function(id) {
            if (grepl("^hmdb:", id, ignore.case = TRUE)) {
                paste0("hmdb:", sub("^[^:]+:", "", id))
            } else if (grepl("^CHEBI:", id, ignore.case = TRUE)) {
                paste0("chebi:", sub("^[^:]+:", "", id))
            } else if (grepl("^LIPIDMAPS:", id, ignore.case = TRUE)) {
                paste0("LIPIDMAPS:", sub("^[^:]+:", "", id))
            } else if (grepl("^HMDB", id, ignore.case = TRUE)) {
                paste0("hmdb:", id)
            } else if (grepl("^LM", id)) {
                paste0("LIPIDMAPS:", id)
            } else {
                id
            }
        }, character(1))
        paste(unique(ids), collapse = "; ")
    }

    data.frame(
        ID = matched$ID,
        Match = TRUE,
        observed_mz = matched$observed_mz,
        Reference_mz = matched$expected_mz,
        Error = matched$ppm_error,
        Adduct = matched$adduct,
        Formula = matched$formula,
        FormulaStatus = matched$formula_status,
        Exactmass = matched$neutral_mass,
        Isomers = matched$metabolite_ids,
        InchiKeys = matched$inchikeys,
        IsomerNames = matched$metabolite_names,
        Isomers_IDs = vapply(
            matched$metabolite_ids,
            standardise_ids,
            character(1)
        ),
        Ramp_IDs = matched$ramp_ids,
        Score = matched$score,
        MassScore = matched$mass_score,
        ChemicalScore = matched$chemical_score,
        WaterLoss = matched$water_loss,
        DehydrationSites = matched$dehydration_sites,
        MaxWaterLosses = matched$max_water_losses,
        DehydrationStatus = matched$dehydration_status,
        StructureScore = matched$structure_score,
        StructureStatus = matched$structure_status,
        StructureRuleEvidence = matched$structure_rule_evidence,
        StructureEvidence = matched$structure_evidence,
        PositiveSiteAtoms = matched$positive_site_atoms,
        NegativeSiteAtoms = matched$negative_site_atoms,
        AlkaliSiteAtoms = matched$alkali_site_atoms,
        ReactiveSiteScore = matched$reactive_site_score,
        ReactiveSiteStatus = matched$reactive_site_status,
        IsotopeScore = matched$isotope_score,
        AdductNetworkScore = matched$adduct_network_score,
        MALDIMatrix = matched$maldi_matrix,
        RuleClass = matched$rule_class,
        RuleSource = matched$rule_source,
        ReactiveGroup = matched$reactive_group,
        MinimumReactiveSites = matched$min_reactive_sites,
        MatrixTargetGroups = matched$matrix_target_groups,
        MatrixTargetSites = matched$matrix_target_sites,
        MatrixTargetStatus = matched$matrix_target_status,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}


.annotation_terms <- function(x) {
    terms <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE))
    unique(terms[!is.na(terms) & nzchar(terms) & terms != "No Annotation"])
}

#' Refines and reduces m/z annotations
#'
#' Used to subset dataset to only include annotations that have n number of
#' entries.
#' For example some peaks can have multiple annotations. Peaks which have above
#' n number of annotations assigned will be removed from Bioconductor
#' experiment.
#'
#' @param obj A SummarizedExperiment, including SingleCellExperiment or
#' SpatialExperiment. The selected experiment must have annotations in
#'   `rowData()`.
#' @param assay Primary (`"main"`) or alternative experiment whose `rowData()`
#'   contains the annotations (default = `"main"`).
#' @param n Integer defining the number of entries an annotation can have
#'   assigned. Any higher counts will be removed (default = 1).
#'
#' @return The input container with the selected experiment restricted to
#'   annotated features having at most `n` candidate rows. `n = 0` returns
#'   no features. Stored candidates are counted when available, so semicolons
#'   within one candidate's isomer names do not inflate the count. For legacy
#'   objects without stored candidates, distinct semicolon-separated names
#'   are counted. Unannotated features are excluded. Pixel columns are unchanged.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' head(getRefinedAnnotations(spe, n = 1))
#'
getRefinedAnnotations <- function(obj, assay = "main", n = 1) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) {
        stop("n must be one non-negative integer.", call. = FALSE)
    }
    metadata <- .featureMetadata(obj, assay)
    selected <- .pathwayAnnotationObject(obj, assay)
    store <- .storedData(selected, "mz_annotation")
    if (.annotation_has_current_schema(store$results) ||
        .annotation_has_current_schema(.storedData(selected, "db_3"))) {
        candidates <- .currentAnnotationCandidates(selected)
        counts <- as.integer(table(factor(candidates$feature, levels = rownames(metadata))))
    } else {
        if (!"all_IsomerNames" %in% names(metadata)) {
            stop("No annotation labels or stored candidates were found.", call. = FALSE)
        }
        counts <- lengths(lapply(as.character(metadata$all_IsomerNames), .annotation_terms))
    }
    subset.obj <- subsetMZFeatures(
        data = obj,
        features = rownames(metadata)[counts > 0L & counts <= n],
        assay = assay
    )
    return(subset.obj)
}

########################################################################################################################################################################################################################


#### SpaMTP Find m/z Annotation Functions #####################################################################################################################################################################################

#' Adds in backslashes required to take into account special using grepl such as
#' brackets and +
#'
#' @param inputString Character string requiring additional backslashes
#'
#' @return Character string containing double backslashes around special
#'   features
#' @noRd
add_backslashes_to_specialfeatures <- function(inputString) {
    special <- c("\\", ".", "^", "$", "|", "(", ")", "[", "]", "{", "}", "*", "+", "?")
    vapply(as.character(inputString), function(value) {
        if (is.na(value)) return(NA_character_)
        chars <- strsplit(value, "", fixed = TRUE)[[1]]
        chars[chars %in% special] <- paste0("\\", chars[chars %in% special])
        paste(chars, collapse = "")
    }, character(1), USE.NAMES = FALSE)
}


#' Find Annotation
#'
#' Searches through annotated m/z values to return all which contain the
#' metabolite search term provided
#'
#' @param data Bioconductor experiment containing annotated m/z values.
#' @param metabolite One literal, case-insensitive metabolite search term.
#'   Punctuation is matched as text, not as a regular expression.
#' @param assay Primary or alternative experiment containing the annotated
#'   m/z metadata. Defaults to `"main"`.
#' @param searchExact Boolean value defining if to only return m/z values which
#'   contain the exact match to the metabolite search term (default = FALSE).
#' @param columnName Character string defining the feature-metadata column
#'   where annotations are stored (default = "all_IsomerNames").
#'
#' @return A Data.Frame containing the peak metadata corresponding to the
#'   metabolite search term provided
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' searchAnnotations(spe, "Glucose", searchExact = TRUE)
searchAnnotations <- function(
    data,
    metabolite,
    assay = "main",
    searchExact = FALSE,
    columnName = "all_IsomerNames"
) {
    if (!is.character(metabolite) || length(metabolite) != 1L ||
        is.na(metabolite) || !nzchar(metabolite)) {
        stop("metabolite must be one non-empty search term.", call. = FALSE)
    }
    search_term <- tolower(metabolite)
    feature_metadata <- .featureMetadata(data, assay)
    if (!"mz_names" %in% names(feature_metadata)) {
        feature_metadata$mz_names <- rownames(feature_metadata)
    }

    if (!columnName %in% names(feature_metadata)) {
        compatible_columns <- c(
            "all_IsomerNames", "all_Isomers", "IsomerNames", "annotation"
        )
        available <- compatible_columns[compatible_columns %in% names(feature_metadata)]
        if (identical(columnName, "all_IsomerNames") && length(available)) {
            columnName <- available[[1]]
        } else {
            # Current serialized/mapped objects can retain the versioned annotation
            # store even when an older feature-label column was not propagated.
            annotation_object <- .pathwayAnnotationObject(data, assay)
            annotations <- tryCatch(
                .get_mz_annotations(annotation_object,
                    annotationSource = "current"
                ),
                error = function(e) NULL
            )
            if (!is.null(annotations) &&
                all(c("observed_mz", "IsomerNames") %in% names(annotations))) {
                annotations$mz_name <- .candidateFeatures(
                    annotation_object, annotations
                )
                annotations <- annotations[!is.na(annotations$mz_name), ,
                    drop = FALSE
                ]
                feature_metadata$mz_names <- rownames(feature_metadata)
                annotation_labels <- annotations %>%
                    dplyr::group_by(.data$mz_name) %>%
                    dplyr::summarise(
                        .annotation_search = paste(unique(.data$IsomerNames), collapse = "; "),
                        .groups = "drop"
                    )
                annotation_labels$mz_names <- annotation_labels$mz_name
                feature_metadata <- dplyr::left_join(
                    feature_metadata,
                    annotation_labels[c("mz_names", ".annotation_search")],
                    by = "mz_names"
                )
                rownames(feature_metadata) <- feature_metadata$mz_names
                columnName <- ".annotation_search"
            } else {
                stop(
                    "Annotation column '", columnName, "' is not present in assay '",
                    assay, "', and no current mz_annotation labels are available."
                )
            }
        }
    }

    annotation_values <- as.character(feature_metadata[[columnName]])
    indexs <- which(!is.na(annotation_values) &
        grepl(search_term, tolower(annotation_values), fixed = TRUE))

    if (searchExact) {
        exact <- vapply(annotation_values[indexs], function(annotationList) {
            split_list <- trimws(unlist(strsplit(annotationList, ";")))
            tolower(metabolite) %in% tolower(split_list)
        }, logical(1))
        indexs <- indexs[exact]
    }

    return(feature_metadata[indexs, , drop = FALSE])
}


#' Finds if any metabolite is duplicated across multiple m/z values.
#'
#' @param data Bioconductor experiment containing annotated m/z values.
#' @param assay Primary or alternative experiment containing the annotated
#'   m/z metadata. Defaults to `"main"`.
#'
#' @return Vector of character strings describing metabolites that are assigned
#'   to multiple m/z values
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' findDuplicateAnnotations(spe)
findDuplicateAnnotations <- function(data, assay = "main") {
    all_annotations <- .featureMetadata(data, assay)$all_IsomerNames
    all_terms <- unlist(lapply(as.character(all_annotations), .annotation_terms), use.names = FALSE)
    terms_counts <- table(all_terms)
    return(names(terms_counts[terms_counts > 1]))
}

#' Gets values from a single metadata column for a respective m/z value.
#'
#' @param obj Bioconductor experiment containing annotated m/z values.
#' @param mz Character string specifying the m/z value to return the metadata
#'   for.
#' @param assay Primary or alternative experiment containing the annotated
#'   m/z metadata. Defaults to `"main"`.
#' @param metadataColumn Character string corresponding to the assay feature
#'   metadata column to extract (default = "all_IsomerNames").
#' @param separate Boolean indicating whether to separate the metadata string
#'   and return a vector. Note, if `TRUE` the metadata column should contain
#'   values separate by "; " (default = TRUE).
#'
#' @return Vector of character strings containing the metadata for a specific mz
#'   value.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' getMZMetadata(spe, mz = "ion1")
#  ##### Example for getting metabolite annotation for a m/z value
#' # getMZMetadata(SpaMTP, mz = "mz-100", metadataColumn = "all_IsomerNames")
#'
#' ##### Example for getting metabolite annotation IDs for a m/z value
#' # getMZMetadata(SpaMTP, mz = "mz-100", metadataColumn = "all_Isomers")
getMZMetadata <- function(
    obj,
    mz,
    assay = "main",
    metadataColumn = "all_IsomerNames",
    separate = TRUE
) {
    df <- .featureMetadata(obj, assay)
    if (length(mz) != 1L || is.na(mz)) stop("Supply one feature name in mz.", call. = FALSE)
    if (!metadataColumn %in% names(df)) {
        stop("Unknown feature metadata column: '", metadataColumn, "'.", call. = FALSE)
    }
    ids <- rownames(df)
    position <- match(as.character(mz), ids)
    if (is.na(position) && "mz_names" %in% names(df)) position <- match(mz, df$mz_names)
    if (is.na(position)) stop("Unknown m/z feature: '", mz, "'.", call. = FALSE)
    mz_row <- df[[metadataColumn]][position]
    if (separate) {
        return(trimws(strsplit(as.character(mz_row), split = ";", fixed = TRUE)[[1]]))
    } else {
        return(mz_row)
    }
}

########################################################################################################################################################################################################################

# !! ALL CODE BELOW WAS WRITEN BY Christopher Fitzgerald github https://github.com/ChemCharles !!#


#' Checks if the complete adduct is in the data base, else returns a truncated
#' adduct
#'
#' @param adduct Character string defining the adduct to be checked.
#' @param db DataFrame of the current reference database.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return Character string of the complete or truncated adduct
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # check_and_truncate_adduct_vector(c("M+H"), HMDB_db)
#' @noRd
check_and_truncate_adduct_vector <- function(adduct, db, verbose = TRUE) {
    element_exists <- adduct %in% colnames(db)
    missing_elements <- adduct[!element_exists]
    if (length(missing_elements) > 0) {
        for (missing_element in missing_elements) {
            verbose_message(messageText = paste0(
                "Adduct ",
                missing_element,
                " is not in the DB, it has been removed from the search."
            ), verbose = verbose)
        }
        truncated_adduct <- adduct[element_exists]
        return(truncated_adduct)
    } else {
        return(adduct)
    }
}


#' Filters the provided metabolomic database by polarity and adducts
#'
#' @param adduct Character string defining the adduct to be checked.
#' @param db DataFrame of the current reference database.
#' @param polarity Character string defining the polarity of the adducts
#'   (default = "neg").
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#'
#' @return A filtered reference metabolomic database DataFrame
#'
#' @examples
#' # HMDB_db <- load("data/HMDB_1_names.rds")
#' # db_adduct_filter(HMDB_db,c("M+H"), polarity = "pos")
#' @noRd
db_adduct_filter <- function(db, adduct, polarity = "neg", verbose = TRUE) {
    # only include adducts from either neg or pos polarity
    if (polarity == "neg") {
        neg_adducts_1 <-
            c(
                "M-H2O-H",
                "M-H",
                "M+Na-2H",
                "M+Cl",
                "M+K-2H",
                "M+FA-H",
                "M+Hac-H",
                "M+Br",
                "M+TFA-H",
                "2M-H",
                "2M+FA-H",
                "2M+Hac-H",
                "3M-H"
            )
        pol <- neg_adducts_1
    } else if (polarity == "pos") {
        pos_adducts_1 <-
            c(
                "M+H",
                "M+NH4",
                "M+Na",
                "M+CH3OH+H",
                "M+K",
                "M+ACN+H",
                "M+2Na-H",
                "M+IsoProp+H",
                "M+ACN+Na",
                "M+2K+H",
                "M+DMSO+H",
                "M+2ACN+H",
                "M+IsoProp+Na+H",
                "2M+H",
                "2M+NH4",
                "2M+Na",
                "2M+K",
                "2M+ACN+H",
                "2M+ACN+Na"
            )
        pol <- pos_adducts_1
    } else {
        stop("Invalid polarity. Choose 'neg' or 'pos'.")
    }

    # get rid of spaces in the adduct names
    # in col names
    db <- db %>%
        dplyr::rename_all(~ gsub(" ", "", .))

    # Filter the db by polarity
    db <- db %>%
        dplyr::select(
            tidyr::all_of(c(
                "formula",
                "exactmass",
                "isomers",
                "isomers_inchikey",
                "isomers_names"
            )),
            tidyr::all_of(pol)
        )

    # in adduct entry
    adduct <- gsub(" ", "", adduct)

    adduct <- check_and_truncate_adduct_vector(adduct, db, verbose = verbose)

    db_filtered <- db %>%
        dplyr::select(
            tidyr::all_of(c(
                "formula",
                "exactmass",
                "isomers",
                "isomers_inchikey",
                "isomers_names"
            )),
            tidyr::any_of(adduct)
        ) %>%
        as.data.frame()
    return(db_filtered)
}


#' Checks if a formula contains only the allowed elements
#'
#' @param formula Character string defining t
#' @param allowedElements Vector of character strings defining allowed elements
#'
#' @noRd
is_formula_valid <- function(formula, allowedElements) {
    elements <-
        stringr::str_extract_all(formula, "[A-Z][a-z]*")[[1]] # defines an element as a Uppercase immediately followed by none or more lowercase letters.

    all(elements %in% allowedElements) # Then checks if they are in the vector of allowed elements.
}


#' Filters reference Database to only select natural elements
#'
#' @param df DataFrame of the reference database.
#' @param elements Vector of character strings of elements to include (default =
#'   c("H", "C", "N", "O", "S", "Cl", "Br", "F", "Na", "P", "I", "Si")).
#'
#' @return A refined DataFrame which only includes annotations containing the
#'   specified elements
#'
#' @examples
#' ## Get filtered DB by adduct
#' # db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos")
#'
#' ## Refine to only include natural elements
#' # db_2 <- formula_filter(db_1)
#' @noRd
formula_filter <- function(df, elements = NULL) {
    if (is.null(elements)) {
        elements <- c(
            "H", "C", "N", "O", "S", "Cl", "Br",
            "F", "Na", "P", "I", "Si"
        )
    }

    # Elements to allow
    allowed_elements <- elements

    # Filter rows based on allowed elements
    filtered_df <- df %>%
        dplyr::filter(vapply(
            .data$formula,
            is_formula_valid,
            logical(1),
            allowedElements = allowed_elements
        ))

    return(filtered_df)
}


#' @noRd
calculate_bounds <- function(inputDf) {
    lower_bound <- min(inputDf$mz, na.rm = TRUE)
    upper_bound <- max(inputDf$mz, na.rm = TRUE)

    bounds <-
        list(lower_bound = lower_bound, upper_bound = upper_bound)

    return(bounds)
}


#' @noRd
ppm_error <- function(observedMz, referenceMz, ppm) {
    abs_diff_ppm <-
        abs(observedMz - referenceMz) / abs(referenceMz) * 1e6
    if (abs_diff_ppm <= ppm) {
        return(abs_diff_ppm)
    } else {
        return("Out")
    }
}


#' @noRd
ppm_range_match <- function(observedMz, referenceMz, ppm) {
    abs_diff_ppm <-
        abs(observedMz - referenceMz) / abs(referenceMz) * 1e6
    abs_diff_ppm <= ppm
}


########################################################################################################################################################################################################################

#' Assign custom annotations to m/z values
#'
#' Adds custom metabolite annotations to respective m/z values (ideal for
#' specific matrices such as FMP10).
#' If several features are exactly tied for the nearest mass, all receive the
#' annotation. Features at unequal distances retain the nearest-mass rule.
#' Existing labels are preserved by default, including on unmatched features.
#' Custom labels are display metadata: they do not add scored candidates or
#' change pathway evidence. Use a separate `annotationColumn` to keep manual
#' labels apart from the candidate summaries produced by `annotateSM()`.
#'
#' @param data Bioconductor experiment containing m/z intensity values.
#' @param annotations data.frame containing two columns named 'annotation' and
#'   'mass'. These columns should contain the custom metabolite annotation and
#'   the relative m/z mass respectively.
#' @param assay Primary or alternative experiment whose `rowData()` receives
#'   the annotations (default = "main").
#' @param returnOnlyAnnotated Boolean defining whether to return a
#'   Bioconductor experiment containing only successfully annotated m/z values
#'   (default = FALSE).
#' @param massThreshold Numeric value defining the acceptable threshold
#'   (plus-minus) between the custom annotations and the actual m/z values
#'   contained within the SpaMTP object, in Da (default = 0.05). The boundary
#'   is inclusive; zero keeps exact matches and NULL disables the distance limit.
#' @param annotationColumn Name of the `rowData()` column receiving annotations
#'   (default = "all_IsomerNames").
#' @param append Append distinct labels and metadata values on matched features
#'   (default TRUE). FALSE replaces values only on matched features.
#'
#' @return Bioconductor experiment containing the custom annotations stored in
#'   the feature metadata dataframe. `mass_diff_da` is the absolute mass
#'   difference in Da; `ppm_diff` is that difference divided by reference mass
#'   and multiplied by one million.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' custom <- data.frame(annotation = "candidate A", mass = 500)
#' annotated <- addCustomMZAnnotations(spe, custom)
#' SummarizedExperiment::rowData(annotated)["ion7", ]
addCustomMZAnnotations <- function(
    data,
    annotations,
    assay = "main",
    returnOnlyAnnotated = FALSE,
    massThreshold = 0.05,
    annotationColumn = "all_IsomerNames",
    append = TRUE
) {
    experiment <- .experimentForAssay(data, assay)
    if (!is.data.frame(annotations) || !all(c("annotation", "mass") %in% colnames(annotations)) ||
        !is.numeric(annotations$mass) || any(!is.finite(annotations$mass) | annotations$mass <= 0)) {
        stop(
            "annotations must contain annotation and positive finite numeric mass columns.",
            call. = FALSE
        )
    }
    if (!is.null(massThreshold) && (!is.numeric(massThreshold) || length(massThreshold) != 1L ||
        !is.finite(massThreshold) || massThreshold < 0)) {
        stop(
            "massThreshold must be NULL or a finite non-negative number.",
            call. = FALSE
        )
    }
    if (!is.logical(append) || length(append) != 1L || is.na(append)) {
        stop("append must be TRUE or FALSE.", call. = FALSE)
    }
    if (length(annotationColumn) != 1L || is.na(annotationColumn) ||
        !nzchar(annotationColumn) || annotationColumn %in%
        c("mz", "raw_mz", "mz_names", ".spamtp_feature_id")) {
        stop(
            "annotationColumn must not replace the feature mass or identifier.",
            call. = FALSE
        )
    }
    featureData <- .featureMetadata(data, assay)
    # A custom label does not replace the candidate store. Never recreate
    # missing keys, which could reconnect stale candidates to different rows.
    if (".spamtp_feature_id" %in% names(featureData)) {
        .annotationFeatureIds(featureData)
    }
    masses <- .massValues(featureData, rownames(experiment))
    matches <- .match_custom_masses(annotations$mass, masses, massThreshold)
    annotations <- annotations[matches$reference_row, , drop = FALSE]
    matched <- matches$feature_row
    annotations$true_mzs <- masses[matched]
    annotations$mass_diff_da <- matches$mass_diff_da
    annotations$ppm_diff <- matches$ppm_diff
    annotations[[annotationColumn]] <- as.character(annotations$annotation)
    if (annotationColumn != "annotation") annotations$annotation <- NULL
    groups <- split(seq_len(nrow(matches)), matched)
    featureData$mz_names <- rownames(experiment)
    protected <- c("mz", "raw_mz", "mz_names", ".spamtp_feature_id")
    for (column in setdiff(colnames(annotations), protected)) {
        values <- if (column %in% names(featureData)) {
            as.character(featureData[[column]])
        } else rep("No Annotation", nrow(experiment))
        positions <- as.integer(names(groups))
        additions <- vapply(groups, function(rows) {
            value <- annotations[[column]][rows]
            paste(unique(value[!is.na(value)]), collapse = "; ")
        }, character(1), USE.NAMES = FALSE)
        if (append) {
            additions <- vapply(seq_along(positions), function(i) {
                paste(.annotation_terms(c(values[positions[i]], additions[i])), collapse = "; ")
            }, character(1))
        }
        values[positions] <- additions
        featureData[[column]] <- values
    }
    data <- .setFeatureMetadata(data, featureData, assay)
    if (isTRUE(returnOnlyAnnotated)) {
        selected <- sort(unique(matched))
        data <- subsetMZFeatures(data, rownames(experiment)[selected], assay)
    }
    data
}

# Match references to the nearest measured mass, preserving exact ties and
# numeric row identities. Never reconstruct these links from display strings.
.match_custom_masses <- function(reference, masses, threshold) {
    empty <- data.frame(reference_row = integer(), feature_row = integer(),
        mass_diff_da = numeric(), ppm_diff = numeric())
    valid <- which(is.finite(masses) & masses > 0)
    if (!length(reference) || !length(valid)) return(empty)
    ordered <- sort(unique(masses[valid]))
    features <- split(valid, match(masses[valid], ordered))
    position <- findInterval(reference, ordered)
    matches <- lapply(seq_along(reference), function(i) {
        neighbours <- unique(pmax(1L, pmin(length(ordered), c(position[i], position[i] + 1L))))
        difference <- abs(ordered[neighbours] - reference[i])
        nearest <- neighbours[difference == min(difference)]
        if (!is.null(threshold) && min(difference) > threshold) return(integer())
        sort(as.integer(unlist(features[as.character(nearest)], use.names = FALSE)))
    })
    ref <- rep(seq_along(reference), lengths(matches))
    feature <- as.integer(unlist(matches, use.names = FALSE))
    difference <- abs(reference[ref] - masses[feature])
    data.frame(reference_row = ref, feature_row = feature,
        mass_diff_da = difference, ppm_diff = difference / reference[ref] * 1e6)
}


#' Annotates FMP10 matrix data
#'
#' Adds metabolite annotations to respective m/z values generated using an FMP10
#' matrix. This is done based on a curated FMP10 matrix database. This
#' annotation pass replaces the selected experiment's candidate summaries
#' and stored candidates with the curated matches. Use `addCustomMZAnnotations()`
#' to add display labels without replacing a candidate annotation pass.
#' Exact nearest-mass ties annotate all tied features, with names, scores and
#' stored candidates kept together under the preserved feature identity.
#'
#' @param obj Bioconductor experiment containing m/z intensity values for
#'   annotation. Data should be generated with a FMP10 matrix.
#' @param onlyFmpAdduct Boolean indicating if only metabolites with FMP10+
#'   adducts (`+FMP10`, +`2FMP10`, etc.) should be assigned to m/z values
#'   (default = FALSE).
#' @param addCustomAnnotation data.frame containing addition FMP10 metabolite
#'   annotations that are not in the current FMP10 database. Note: this
#'   data.frame must contain these column c("mass", "annotation", "Adduct",
#'   "Formula", "Isomers", "Isomers_IDs"). If set to NULL, only reference FMP10
#'   database will be used (default = NULL).
#' @param assay Primary or alternative experiment whose `rowData()` receives
#'   the annotations (default = "main").
#' @param returnOnlyAnnotated Boolean defining whether to return a
#'   Bioconductor experiment containing only successfully annotated m/z values
#'   (default = FALSE).
#' @param massThreshold Numeric value defining the acceptable threshold
#'   (plus-minus) between the custom annotations and the actual m/z values
#'   contained within the SpaMTP object, in Da (default = 0.05). The boundary
#'   is inclusive; zero keeps exact matches and NULL disables the distance limit.
#'   Choose this tolerance for the calibrated instrument accuracy; 0.05 Da
#'   can be permissive for high-resolution data. `all_Errors` is reported in ppm.
#' @param annotationColumn Name of the `rowData()` column receiving annotations
#'   (default = "all_IsomerNames").
#' @param database Optional named list of database resources, normally created
#'   by [loadSpaMTPDatabase()].
#' @param databaseVersion SpaMTPdb/RaMP version used for annotation lookup.
#' @param databaseSource Database source; see [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#'
#' @return Bioconductor experiment containing the relative metabolite
#'   annotations stored in the feature metadata dataframe.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' spe <- spamtp_example$experiment
#' # A miniature supplied reference keeps the example offline.
#' reference <- data.frame(
#'     mass = 500, annotation = "example derivative",
#'     Adduct = "M+FMP10", Formula = "C10H10O2", Isomers = "example",
#'     Isomers_IDs = "example:1", IsomerNames = "example derivative"
#' )
#' annotated <- addFMP10Annotations(spe,
#'     database = list(filtered_fmp10 = reference, chem_props = spamtp_example$metabolites)
#' )
#' SummarizedExperiment::rowData(annotated)["ion7", ]
addFMP10Annotations <- function(
    obj, onlyFmpAdduct = FALSE,
    addCustomAnnotation = NULL,
    assay = "main",
    returnOnlyAnnotated = FALSE,
    massThreshold = 0.05,
    annotationColumn = "all_IsomerNames",
    database = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL
) {
    database_resources <- .spamtp_db_bundle(
        c("filtered_fmp10", "chem_props"),
        database = database,
        version = databaseVersion,
        source = match.arg(databaseSource),
        localDir = databaseLocalDir
    )
    filtered_fmp10 <- database_resources$filtered_fmp10
    chem_props <- database_resources$chem_props

    if (!is.null(massThreshold) &&
        (!is.numeric(massThreshold) || length(massThreshold) != 1L ||
            !is.finite(massThreshold) || massThreshold < 0)) {
        stop("massThreshold must be NULL or one non-negative finite number.")
    }

    if (onlyFmpAdduct) {
        filtered_fmp10 <- filtered_fmp10[grepl(x = filtered_fmp10$Adduct, pattern = "FMP10"), ]
    }

    if (!is.null(addCustomAnnotation)) {
        required_columns <- c(
            "mass",
            "annotation",
            "Adduct",
            "Formula",
            "Isomers",
            "Isomers_IDs"
        )

        if (all(required_columns %in% colnames(addCustomAnnotation))) {
            addCustomAnnotation$IsomerNames <- addCustomAnnotation$annotation
            addCustomAnnotation <- addCustomAnnotation[c(required_columns, "IsomerNames")]
            filtered_fmp10 <- dplyr::bind_rows(filtered_fmp10, addCustomAnnotation)
        } else {
            stop("data.frame provided in `addCustomAnnotation` does not contain the correct column names. Column names must be: c('mass', 'annotation', 'Adduct', 'Formula', 'Isomers', 'Isomers_IDs')")
        }
    }


    featureMetadata <- .featureMetadata(obj, assay)
    featureMetadata$.spamtp_feature_id <- .annotationFeatureIds(
        featureMetadata,
        initialize = TRUE
    )
    if (!is.character(annotationColumn) || length(annotationColumn) != 1L ||
        is.na(annotationColumn) || !nzchar(annotationColumn) ||
        annotationColumn %in% c("mz", "raw_mz", "mz_names", ".spamtp_feature_id")) {
        stop("annotationColumn must not replace the feature mass or identifier.", call. = FALSE)
    }
    if (!is.numeric(filtered_fmp10$mass) ||
        any(!is.finite(filtered_fmp10$mass) | filtered_fmp10$mass <= 0)) {
        stop("FMP10 reference masses must be positive finite numbers.", call. = FALSE)
    }
    masses <- .massValues(featureMetadata, rownames(featureMetadata))
    matches <- .match_custom_masses(filtered_fmp10$mass, masses, massThreshold)
    db_3 <- filtered_fmp10[matches$reference_row, , drop = FALSE]
    db_3$feature <- rownames(featureMetadata)[matches$feature_row]
    db_3$feature_id <- featureMetadata$.spamtp_feature_id[matches$feature_row]
    db_3$observed_mz <- masses[matches$feature_row]

    # Preserve the historical table, but also map source identifiers onto the
    # current RaMP compound graph and store scored, versioned candidates for all
    # new pathway/network functions.
    alternatives <- if (inherits(obj, "SingleCellExperiment")) {
        SingleCellExperiment::altExpNames(obj)
    } else character()
    stores <- .storedData(obj, "modality_annotations") %||% list()
    previous <- .storedData(obj, "mz_annotation")
    previous_legacy <- .storedData(obj, "db_3")
    previous_assay <- previous$metadata$assay %||% "main"
    previous_key <- if (previous_assay %in% alternatives) {
        previous_assay
    } else "main"
    if (!is.null(previous) || !is.null(previous_legacy)) {
        stores[[previous_key]] <- list(
            mz_annotation = previous, db_3 = previous_legacy
        )
    }
    obj <- .setStoredData(obj, "db_3", db_3)
    current_annotations <- .fmp10_current_annotations(
        db_3,
        massThreshold = massThreshold,
        chemicalProperties = chem_props
    )
    obj <- .setStoredData(obj, "mz_annotation", .annotation_store(
        current_annotations,
        metadata = list(
            engine = "curated-fmp10-v2",
            assay = assay,
            mass_tolerance_da = massThreshold,
            ramp_version = .annotation_ramp_version(chem_props),
            database = paste0(
                "curated FMP10 panel mapped to ", .spamtp_db_label(chem_props)
            ),
            mapped_candidates = sum(nzchar(current_annotations$Ramp_IDs))
        )
    ))
    key <- if (assay %in% alternatives) assay else "main"
    stores[[key]] <- list(
        mz_annotation = .storedData(obj, "mz_annotation"), db_3 = db_3
    )
    obj <- .setStoredData(obj, "modality_annotations", stores)

    annotation_summary <- current_annotations %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::summarise(
            all_IsomerNames = paste(.data$IsomerNames, collapse = "; "),
            all_Isomers = paste(.data$Isomers, collapse = "; "),
            all_Isomers_IDs = paste(.data$Isomers_IDs, collapse = "; "),
            all_Adducts = paste(.data$Adduct, collapse = "; "),
            all_Formulas = paste(.data$Formula, collapse = "; "),
            all_Errors = paste(round(.data$Error, 4), collapse = "; "),
            all_Ramp_IDs = paste(
                unique(.data$Ramp_IDs[nzchar(.data$Ramp_IDs)]),
                collapse = "; "
            ),
            all_Scores = paste(round(.data$Score, 4), collapse = "; "),
            .groups = "drop"
        )
    matched <- match(
        featureMetadata$.spamtp_feature_id,
        annotation_summary$feature_id
    )
    for (column in setdiff(names(annotation_summary), "feature_id")) {
        value <- annotation_summary[[column]][matched]
        value[is.na(value)] <- "No Annotation"
        featureMetadata[[column]] <- value
    }

    featureMetadata[[annotationColumn]] <- featureMetadata$all_IsomerNames
    featureMetadata$mz_names <- rownames(featureMetadata)
    obj <- .setFeatureMetadata(obj, featureMetadata, assay)
    if (isTRUE(returnOnlyAnnotated)) {
        obj <- subsetMZFeatures(obj, rownames(featureMetadata)[!is.na(matched)], assay)
    }
    return(obj)
}


.fmp10_current_annotations <- function(
    db3, massThreshold,
    chemicalProperties
) {
    required <- c(
        "mass", "observed_mz", "Adduct", "Formula", "Isomers",
        "Isomers_IDs", "IsomerNames"
    )
    missing <- setdiff(required, names(db3))
    if (length(missing)) {
        stop(
            "FMP10 annotation result is missing required column(s): ",
            paste(missing, collapse = ", ")
        )
    }
    if (!all(c("ramp_id", "chem_source_id") %in% names(chemicalProperties))) {
        stop("chemicalProperties must contain ramp_id and chem_source_id columns.")
    }

    result <- tidyr::separate_rows(
        as.data.frame(db3, stringsAsFactors = FALSE),
        "Isomers_IDs",
        sep = "\\s*;\\s*"
    )
    result$.source_key <- toupper(trimws(as.character(result$Isomers_IDs)))

    source_map <- unique(data.frame(
        .source_key = toupper(trimws(as.character(chemicalProperties$chem_source_id))),
        ramp_id = toupper(trimws(as.character(chemicalProperties$ramp_id))),
        stringsAsFactors = FALSE
    ))
    source_map <- source_map[
        nzchar(source_map$.source_key) & grepl("^RAMP_C_", source_map$ramp_id), ,
        drop = FALSE
    ]
    ramp_lookup <- vapply(
        split(source_map$ramp_id, source_map$.source_key),
        function(x) paste(unique(x), collapse = "; "),
        character(1)
    )
    result$Ramp_IDs <- unname(ramp_lookup[result$.source_key])
    result$Ramp_IDs[is.na(result$Ramp_IDs)] <- ""

    result$observed_mz <- suppressWarnings(as.numeric(result$observed_mz))
    result$Reference_mz <- suppressWarnings(as.numeric(result$mass))
    absolute_error <- abs(result$observed_mz - result$Reference_mz)
    result$Error <- absolute_error / result$Reference_mz * 1e6
    if (is.null(massThreshold) || massThreshold == 0) {
        result$MassScore <- rep(1, nrow(result))
    } else {
        result$MassScore <- exp(-0.5 * (absolute_error / (massThreshold / 3))^2)
    }
    result$ChemicalScore <- rep(1, nrow(result))
    result$IsotopeScore <- rep(NA_real_, nrow(result))
    result$AdductNetworkScore <- rep(NA_real_, nrow(result))
    result$Score <- result$MassScore
    result$Match <- rep(TRUE, nrow(result))
    feature <- if ("feature_id" %in% names(result)) {
        result$feature_id
    } else {
        result$observed_mz
    }
    result$ID <- match(feature, unique(feature))
    # For derivatised FMP10 entries this is the curated panel mass, not a
    # re-derived neutral mass; retain it as the best available exact reference.
    result$Exactmass <- result$Reference_mz
    if (!"InchiKeys" %in% names(result)) result$InchiKeys <- rep(NA_character_, nrow(result))
    result$.source_key <- NULL

    preferred <- c(
        "ID", "Match", "observed_mz", "Reference_mz", "Error", "Adduct",
        "Formula", "Exactmass", "Isomers", "InchiKeys", "IsomerNames",
        "Isomers_IDs", "Ramp_IDs", "Score", "MassScore", "ChemicalScore",
        "IsotopeScore", "AdductNetworkScore"
    )
    result[, c(preferred, setdiff(names(result), preferred)), drop = FALSE]
}


#' Annotates vector of m/z values
#'
#' This function assigns each valid m/z peak with one/multiple metabolite names
#' based on the mass difference between the observed value and the theoretical
#' value documented in the reference database.
#' This function is to be used when dealing with large datasets as a
#' preprocessing step. Users can annotate m/z values first and then subset their
#' data accordinly before loading it into a Bioconductor experiment.
#'
#' @param mzs Vector containing m/z values for annotation.
#' @param db Reference metabolite dataset in the form of a data.frame. SpaMTP
#'   loads versioned reference resources through SpaMTPdb when `db = NULL`.
#'   Supply a data.frame for a custom reference, or a pre-built `index`
#'   through `...`.
#' @param ppm Mass tolerance in ppm (default 5). Set to zero
#'   for exact numerical matches.
#' @param adducts Optional adduct names/notations; see `adductRules()`. If
#'   `NULL`, use the complete rule space selected from `maldiMatrix`, or all
#'   validated general rules for the selected polarity when no matrix is given.
#' @param polarity Character string defining the ion mode. When `NULL`, infer
#'   it from a supplied index or MALDI matrix profile, otherwise use positive.
#' @param verbose Boolean indicating whether to show the message. If TRUE the
#'   message will be show, else the message will be suppressed (default = TRUE).
#' @param maldiMatrix Optional MALDI matrix/reagent profile used for automatic
#'   rule selection. The `adducts` argument remains an optional restriction.
#' @param databaseVersion SpaMTPdb/RaMP version used when `db = NULL`.
#' @param databaseSource Database source used when `db = NULL`; see
#'   [loadSpaMTPDatabase()].
#' @param databaseLocalDir Optional staged SpaMTPdb resource directory.
#' @param ... Additional indexed annotation/scoring arguments passed to
#'   `annotateTable()`, such as `index`, `rules`, or `ms1_spectrum`.
#'
#' @returns A data.frame containing all successfully annotated m/z values, with
#'   their corresponding annotation.
#' @export
#'
#' @examples
#' data(spamtp_example)
#' annotateBigData(c(149.04445, 181.07066),
#'     db = spamtp_example$metabolites,
#'     adducts = "M+H", ppm = 5, verbose = FALSE
#' )
annotateBigData <- function(
    mzs, db = NULL, ppm = 5, adducts = NULL,
    polarity = NULL,
    verbose = TRUE, maldiMatrix = NULL,
    databaseVersion = "latest",
    databaseSource = c("auto", "spamtpdb"),
    databaseLocalDir = NULL, ...
) {
    mz_df <- data.frame(mz = unique(mzs))
    mz_df$row_id <- seq_along(mz_df[["mz"]])
    mz_df <- mz_df[c("row_id", "mz")]

    annotations <- annotateTable(
        mz_df,
        db = db, ppm = ppm, adducts = adducts,
        polarity = polarity, verbose = verbose,
        maldiMatrix = maldiMatrix, databaseVersion = databaseVersion,
        databaseSource = match.arg(databaseSource),
        databaseLocalDir = databaseLocalDir, ...
    )

    result_df <- annotations %>%
        dplyr::group_by(.data$observed_mz) %>%
        dplyr::summarise(
            all_IsomerNames = paste(.data$IsomerNames, collapse = "; "),
            all_Isomers = paste(.data$Isomers, collapse = "; "),
            all_Isomers_IDs = paste(.data$Isomers_IDs, collapse = "; "),
            all_Adducts = paste(unique(.data$Adduct), collapse = "; "),
            all_Formulas = paste(unique(.data$Formula), collapse = "; "),
            all_Errors = paste(round(.data$Error, 4), collapse = "; "),
            all_Scores = paste(round(.data$Score, 4), collapse = "; "),
            all_Ramp_IDs = paste(
                unique(.data$Ramp_IDs[!is.na(.data$Ramp_IDs)]),
                collapse = "; "
            )
        )

    result_df <- data.frame(result_df)
    result_df$mz_names <- if (nrow(result_df)) {
        paste0("mz-", result_df$observed_mz)
    } else {
        character()
    }
    result_df$present <- rep(TRUE, nrow(result_df))
    rownames(result_df) <- seq_len(nrow(result_df))
    return(result_df)
}
