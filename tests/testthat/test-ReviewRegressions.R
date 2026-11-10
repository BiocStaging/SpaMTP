# Development assistance: see inst/COPYRIGHTS.

test_that("spatial plotting accepts unnamed pixels without changing the input", {
    object <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = rbind(`mz-100` = 1:6, `mz-200` = 6:1)),
        spatialCoords = cbind(x = rep(1:3, 2), y = rep(1:2, each = 3))
    )
    before <- object
    plots <- plotSpatialFeature(object, c("mz-100", "mz-200"), combine = FALSE)
    expect_equal(plots[[1]]$data$intensity, as.numeric(1:6))
    expect_equal(plots[[2]]$data$intensity, as.numeric(6:1))
    expect_equal(plots[[1]]$data$x, rep(1:3, 2))
    expect_equal(plots[[1]]$data$y, rep(1:2, each = 3))
    expect_s3_class(plotSpatialFeature(object, c("mz-100", "mz-200")), "ggplot")
    expect_identical(object, before)
    expect_error(.nativeCoordinates(object), "unique, non-missing names")
})

test_that("PCA stores a variance plot without opening a graphics device", {
    object <- nativeFixture()
    device <- grDevices::dev.cur()
    result <- runMetabolicPCA(object,
        npcs = 2, show_variance_plot = TRUE,
        verbose = FALSE
    )
    expect_identical(grDevices::dev.cur(), device)
    expect_s3_class(S4Vectors::metadata(result)$reduction_plots$pca, "ggplot")
})

test_that("regional enrichment stays silent and preserves explicit parallel parameters", {
    ids <- paste0("RAMP_G_", seq_len(40))
    db <- list(
        source_df = data.frame(rampId = ids, commonName = ids),
        analytehaspathway = data.frame(rampId = ids[1:8], pathwayRampId = "P1"),
        pathway = data.frame(pathwayRampId = "P1", pathwayName = "one")
    )
    index <- buildPathwayIndex(db, gene_mapping = "ramp", organism = "custom")
    ranks <- withr::with_seed(11, matrix(stats::rnorm(40),
        dimnames = list(ids, "A")
    ))
    expect_output(result <- findRegionalPathways(
        ranks = list(genes = ranks),
        pathway_index = index, min_path_size = 3, verbose = FALSE
    ), NA)
    expect_equal(nrow(result), 1L)
    parameter <- BiocParallel::SerialParam(progressbar = FALSE)
    expect_identical(findRegionalPathways(
        ranks = list(genes = ranks),
        pathway_index = index, min_path_size = 3, verbose = FALSE,
        BPPARAM = parameter
    ), result)
})

test_that("METASPACE downloads respect BPPARAM and leave session options alone", {
    options_before <- options()
    parameter <- BiocParallel::SerialParam(progressbar = FALSE)
    testthat::local_mocked_bindings(metaspaceClient = function(...) {
        list(
            get_results = function(...) {
                data.frame(
                    mz = c(200, 100), formula = c("B", "A"),
                    adduct = "M+H", chemMod = "", neutralLoss = "",
                    isotopeImages = I(list(list(list(url = "B")), list(list(url = "A"))))
                )
            },
            get_ion_image_pixels = function(sum_formula, ...) {
                stop("Image metadata should not be queried again")
            }
        )
    }, .metaspace_image_pixels = function(img_info, ...) {
        matrix(if (img_info$url == "A") 1 else 2, 2, 2)
    })
    result <- getMetaspace("offline-test", verbose = FALSE, BPPARAM = parameter)
    expect_identical(names(result$images), c("100.00000", "200.00000"))
    expect_equal(result$images[[1]], matrix(1, 2, 2))
    expect_identical(options(), options_before)
})

test_that("formula helpers use complete element counts and raster export preserves colours", {
    expect_equal(.element_count("CH3OH", "H"), 4L)
    expect_equal(.element_count("NaCl", "C"), 0L)
    expect_equal(
        matrixIonTable("DHB")$mz[
            matrixIonTable("DHB")$ion == "[M]+."
        ],
        unname(MetaboCoreUtils::calculateMass("C7H6O4")) - .spamtp_mass_constants()$electron
    )
    path <- tempfile(fileext = ".png")
    on.exit(unlink(path))
    input <- as.raster(matrix(c("red", "blue", "green", "#00000000", "white", "black"), 2))
    .writeRasterPNG(input, path)
    result <- as.raster(png::readPNG(path))
    expect_identical(
        grDevices::col2rgb(as.vector(result), alpha = TRUE),
        grDevices::col2rgb(as.vector(input), alpha = TRUE)
    )
})
