test_that("visualisePathways returns an informative plot for an empty min_n result", {
    pathway_df <- data.frame(
        analytes_in_pathways = c(1L, 2L),
        p_val = c(0.01, 0.02)
    )

    plot <- visualisePathways(
        object = NULL,
        pathwayDf = pathway_df,
        minN = 3,
        verbose = FALSE
    )

    expect_s3_class(plot, "ggplot")
    expect_identical(plot$labels$title, "No pathways to visualise")
})

test_that("visualisePathways returns an informative plot after p-value filtering", {
    pathway_df <- data.frame(
        pathway_name = "Example pathway",
        pathway_id = "example:1",
        p_val = 0.2,
        analytes_in_pathways = 3L,
        total_in_pathways = 10L,
        adduct_info = ""
    )

    plot <- visualisePathways(
        object = NULL,
        pathwayDf = pathway_df,
        minN = 3,
        pValThreshold = 0.1,
        verbose = FALSE
    )

    expect_s3_class(plot, "ggplot")
    expect_identical(plot$labels$title, "No pathways to visualise")
})

test_that("pathway displays retain supplied statistics and wrap long names", {
    pathway_df <- data.frame(
        pathway_name = "Oligodendrocyte specification and differentiation",
        pathway_id = "WP4304", p_val = 0.000123, fdr = 0.019,
        analytes_in_pathways = 5L, total_in_pathways = 8L
    )
    # Gene results and pathwayAllInfo = FALSE need no ion-image information.
    plot <- visualisePathways(NULL, pathway_df, verbose = FALSE, labelWidth = 20)
    expect_no_warning(ggplot2::ggplotGrob(plot))
    labels <- plot$scales$get_scales("x")$labels
    expect_true(endsWith(labels, "\nWP4304"))
    expect_gt(length(strsplit(labels, "\n")[[1]]), 2L)
    expect_equal(plot$data$p_val, pathway_df$p_val)
    expect_equal(plot$data$fdr, pathway_df$fdr)
    expect_match(plot$data$.statistics, "P = 0.00012\nFDR = 0.019", fixed = TRUE)
    expect_equal(ggplot2::ggplot_build(plot)$data[[1]]$y,
        pathway_df$total_in_pathways)
    expect_identical(plot$theme$legend.position, "bottom")

    pathway_df$fdr <- NULL
    plot <- visualisePathways(NULL, pathway_df, verbose = FALSE, labelWidth = NULL)
    expect_false(grepl("FDR", plot$data$.statistics))
    expect_identical(plot$scales$get_scales("x")$labels,
        paste0(pathway_df$pathway_name, "\n", pathway_df$pathway_id))
    expect_error(visualisePathways(NULL, pathway_df, labelWidth = 0), "labelWidth")
    expect_error(visualisePathways(NULL, pathway_df, textSize = NA_real_), "textSize")
})

test_that("identical pathway memberships render a finite dendrogram", {
    pathways <- data.frame(
        pathway_name = c("First pathway", "Second pathway"),
        pathway_id = c("S1", "S2"), p_val = .01,
        analytes_in_pathways = 1L, total_in_pathways = 1L
    )
    db <- list(
        pathway = data.frame(pathwayRampId = c("P1", "P2"),
            sourceId = c("S1", "S2")),
        analytehaspathway = data.frame(pathwayRampId = c("P1", "P2"),
            rampId = c("A", "A"))
    )
    expect_no_warning({
        plot <- visualisePathways(NULL, pathways, minN = 1,
            database = db, verbose = FALSE)
        ggplot2::ggplotGrob(plot)
    })
})
