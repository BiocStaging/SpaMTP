# Development assistance: see inst/COPYRIGHTS.

pathwayFeatureIdentityFixture <- function(equal_mass = FALSE) {
    mass <- if (equal_mass) c(100.1, 100.1) else c(100.1, 200.2)
    labels <- if (equal_mass) c("ion_a", "ion_b") else {
        c("mz-100.1000", "mz-200.2000")
    }
    values <- matrix(seq_len(8), nrow = 2,
        dimnames = list(labels, paste0("pixel", seq_len(4)))
    )
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = values),
        rowData = S4Vectors::DataFrame(
            mz = mass, .spamtp_feature_id = c("first", "second")
        ),
        spatialCoords = cbind(x = seq_len(4), y = 0)
    )
    ids <- c("RAMP_C_1", "RAMP_C_2")
    candidates <- data.frame(
        observed_mz = mass, feature_id = c("first", "second"),
        mz_name = labels, Adduct = "M+H", Ramp_IDs = ids,
        Score = 1, MassScore = 1, ChemicalScore = 1,
        IsotopeScore = NA_real_, AdductNetworkScore = NA_real_, Error = 0
    )
    S4Vectors::metadata(object)$mz_annotation <- list(
        metadata = list(
            schema_version = 2L, engine = "indexed-chemical-v2",
            ramp_version = "example", assay = "main", raw_mz_column = "mz"
        ),
        results = candidates
    )
    chemical <- data.frame(
        ramp_id = ids, chem_source_id = c("example:1", "example:2"),
        common_name = c("compound A", "compound B")
    )
    database <- list(chem_props = chemical, source_df = data.frame(
        rampId = ids, sourceId = chemical$chem_source_id,
        commonName = chemical$common_name
    ))
    list(object = object, database = database, ids = ids)
}

test_that("compound assays retain padded, renamed and subset features", {
    fixture <- pathwayFeatureIdentityFixture()
    object <- fixture$object
    original <- createPathwayAssay(object, "metabolites",
        assayName = "counts", database = fixture$database,
        outputExperiment = "compounds", verbose = FALSE
    )
    expected <- SummarizedExperiment::assay(object, "counts")
    rownames(expected) <- fixture$ids
    actual <- SingleCellExperiment::altExp(original, "compounds")
    expect_equal(as.matrix(SummarizedExperiment::assay(actual)), expected)

    renamed <- object[2:1, ]
    rownames(renamed) <- c("renamed_second", "renamed_first")
    expect_identical(
        .annotationStatisticsCandidates(renamed, "main",
            fixture$database$source_df
        ),
        list(renamed_first = "RAMP_C_1", renamed_second = "RAMP_C_2")
    )
    mapped <- createPathwayAssay(renamed, "metabolites",
        assayName = "counts", database = fixture$database,
        outputExperiment = "compounds", verbose = FALSE
    )
    compounds <- SingleCellExperiment::altExp(mapped, "compounds")
    expect_equal(as.matrix(SummarizedExperiment::assay(compounds)), expected)
    expect_identical(
        S4Vectors::metadata(compounds)$pathway_mapping$original_features,
        stats::setNames(list("renamed_first", "renamed_second"), fixture$ids)
    )
    retained <- createPathwayAssay(renamed[1, ], "metabolites",
        assayName = "counts", database = fixture$database,
        outputExperiment = "compounds", verbose = FALSE
    )
    expect_equal(as.matrix(SummarizedExperiment::assay(
        SingleCellExperiment::altExp(retained, "compounds")
    )), expected[2, , drop = FALSE])
})

test_that("compound identity follows equal-mass rows in an alternative assay", {
    fixture <- pathwayFeatureIdentityFixture(equal_mass = TRUE)
    object <- fixture$object
    reordered <- object[2:1, ]
    rownames(reordered) <- c("second ion", "first ion")
    SingleCellExperiment::altExp(object, "ions") <- reordered
    mapped <- createPathwayAssay(object, "metabolites", assay = "ions",
        assayName = "counts", database = fixture$database,
        outputExperiment = "compounds", verbose = FALSE
    )
    expected <- SummarizedExperiment::assay(fixture$object, "counts")
    rownames(expected) <- fixture$ids
    compounds <- SingleCellExperiment::altExp(mapped, "compounds")
    expect_equal(as.matrix(SummarizedExperiment::assay(compounds)), expected)
    expect_identical(
        S4Vectors::metadata(compounds)$pathway_mapping$original_features,
        stats::setNames(list("first ion", "second ion"), fixture$ids)
    )
})

test_that("older compound stores use unique masses and reject ambiguous ones", {
    for (equal_mass in c(FALSE, TRUE)) {
        fixture <- pathwayFeatureIdentityFixture(equal_mass)
        object <- fixture$object
        SummarizedExperiment::rowData(object)$.spamtp_feature_id <- NULL
        store <- S4Vectors::metadata(object)$mz_annotation
        store$results$feature_id <- NULL
        S4Vectors::metadata(object)$mz_annotation <- store
        rownames(object) <- c("new_first", "new_second")
        if (equal_mass) {
            expect_error(createPathwayAssay(object, "metabolites",
                assayName = "counts", database = fixture$database,
                verbose = FALSE
            ), "ambiguous feature masses")
        } else {
            expect_silent(mapped <- createPathwayAssay(object, "metabolites",
                assayName = "counts", database = fixture$database,
                outputExperiment = "compounds", verbose = FALSE
            ))
            expect_identical(rownames(SingleCellExperiment::altExp(
                mapped, "compounds"
            )), fixture$ids)
        }
    }
})

test_that("regional and network joins retain renamed metabolite features", {
    fixture <- pathwayFeatureIdentityFixture()
    object <- fixture$object
    object$region <- rep(c("A", "B"), each = 2)
    database <- fixture$database
    database$analyte <- data.frame(rampId = fixture$ids)
    database$analytehaspathway <- data.frame(
        rampId = fixture$ids, pathwayRampId = c("P1", "P2")
    )
    database$pathway <- data.frame(
        pathwayRampId = c("P1", "P2"),
        pathwayName = c("First pathway", "Second pathway"),
        sourceId = c("WP1", "WP2"), type = "WikiPathways"
    )
    database$ramp_db_metadata <- list(ramp_version = "example")
    database$ramp_wikipathway <- list(list(
        id = "WP1", title = "First pathway",
        metabolEdges = data.frame(
            src = fixture$ids[1], dest = fixture$ids[2],
            directed = 1L, reaction_type = 1L
        )
    ))
    index <- buildPathwayIndex(database,
        gene_mapping = "ramp", organism = "custom"
    )
    de <- data.frame(
        gene = rownames(object), cluster = "A", avg_log2FC = c(2, -1),
        p_val_adj = c(0.01, 0.02)
    )
    original <- findRegionalPathways(object, ident = "region",
        DE.list = list(metabolites = de), analyte_types = "metabolites",
        min_path_size = 1, max_path_size = 1, database = database,
        pathway_index = index, organism = "custom", verbose = FALSE
    )
    rownames(object) <- c("renamed_first", "renamed_second")
    de$gene <- rownames(object)
    renamed <- findRegionalPathways(object, ident = "region",
        DE.list = list(metabolites = de), analyte_types = "metabolites",
        min_path_size = 1, max_path_size = 1, database = database,
        pathway_index = index, organism = "custom", verbose = FALSE
    )
    expect_equal(nrow(renamed), 2L)
    expect_equal(renamed, original)
    regional <- renamed[renamed$pathwayRampId == "P1", , drop = FALSE]
    path <- pathwayNetworkPlots(object, ident = "region",
        regpathway = regional, DE.list = list(metabolites = de),
        analyte_types = "metabolites", database = database,
        pathway_index = index, organism = "custom",
        path = withr::local_tempdir(), verbose = FALSE
    )
    expect_true(file.exists(path))
    expect_gte(attr(path, "network_summary")$measured_compounds, 1L)
    html <- paste(readLines(path, warn = FALSE), collapse = "\n")
    payload <- sub(".*const payload = ([^\n]+);\n.*", "\\1", html)
    wire <- jsonlite::fromJSON(payload, simplifyVector = FALSE)
    nodes <- wire$networks[[1]][[1]]$nodes
    first <- Filter(function(node) node$id == fixture$ids[1], nodes)
    expect_length(first, 1L)
    expect_true(first[[1]]$annotated)
})

test_that("annotation search preserves equal-mass keys and modality selection", {
    fixture <- pathwayFeatureIdentityFixture(equal_mass = TRUE)
    object <- fixture$object
    SummarizedExperiment::rowData(object)$mz_names <- rownames(object)
    store <- S4Vectors::metadata(object)$mz_annotation
    store$results$IsomerNames <- c("candidate A", "candidate B")
    S4Vectors::metadata(object)$mz_annotation <- store
    reordered <- object[2:1, ]
    rownames(reordered) <- c("second ion", "first ion")
    expect_identical(rownames(searchAnnotations(reordered, "candidate A")),
        "first ion"
    )
    expect_identical(rownames(searchAnnotations(reordered, "candidate B")),
        "second ion"
    )
    expect_equal(nrow(searchAnnotations(reordered[1, ], "candidate A")), 0L)
    SingleCellExperiment::altExp(object, "ions") <- reordered
    expect_identical(rownames(searchAnnotations(object, "candidate B",
        assay = "ions"
    )), "second ion")
    store$results$feature_id <- NULL
    S4Vectors::metadata(object)$mz_annotation <- store
    SummarizedExperiment::rowData(object)$.spamtp_feature_id <- NULL
    expect_error(searchAnnotations(object, "candidate A"),
        "ambiguous feature masses"
    )
})
