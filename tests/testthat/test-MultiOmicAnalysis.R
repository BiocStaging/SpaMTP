test_that("multiOmicIntegration restores future connection diagnostics", {
    previous <- getOption("future.connections.onMisuse")

    suppressWarnings(
        try(
            multiOmicIntegration(
                NULL,
                reductionList = list("missing-a", "missing-b"),
                dimsList = list(1, 1)
            ),
            silent = TRUE
        )
    )

    expect_identical(getOption("future.connections.onMisuse"), previous)
})
