metaspace_review_client <- function() {
    client <- metaspaceClient()
    annotation <- list(sumFormula = "C6H12O6", adduct = "+H", mz = 181.070664,
        msmScore = 0.9, fdrLevel = 0.1, rhoSpatial = 0.9, rhoSpectral = 0.9,
        rhoChaos = 0.9, offSample = FALSE, possibleCompounds = list(list(name = "hexose")),
        isotopeImages = list(
            list(mz = 181.070664, url = "https://example.invalid/image.png", minIntensity = 0, maxIntensity = 100),
            list(mz = 182.074019, url = NULL, minIntensity = 0, maxIntensity = 10)))
    second <- annotation
    second$mz <- 181.070665
    second$isotopeImages[[1]]$url <- NULL
    second$rhoSpectral <- NULL
    environment(client$get_results)$query_graphql <- function(query, variables = list()) {
        if (grepl("allMolecularDBs", query, fixed = TRUE)) {
            return(list(allMolecularDBs = list(list(id = "test-db", name = "HMDB", version = "v4"))))
        }
        list(allAnnotations = if (is.null(variables$offset) || variables$offset == 0) list(annotation, second) else list())
    }
    client
}

test_that("METASPACE null URLs remain list entries through pagination and image lookup", {
    client <- metaspace_review_client()
    ann <- client$get_results("fixture")
    expect_equal(nrow(ann), 2L)
    expect_identical(names(ann), c("formula", "adduct", "chemMod", "neutralLoss", "mz", "msm", "fdr", "rhoSpatial", "rhoSpectral", "moc", "offSample", "intensity", "moleculeNames", "isotopeImages"))
    expect_true(is.list(ann$isotopeImages))
    expect_length(ann$isotopeImages[[1]], 2L)
    expect_null(ann$isotopeImages[[1]][[2]]$url)
    expect_true(is.na(ann$rhoSpectral[2]))
    local_mocked_bindings(.metaspace_image_pixels = function(img_info, ...) {
        if (is.null(img_info$url)) NULL else matrix(1:6, 2)
    })
    imgs <- client$all_annotation_images("fixture", only_first_isotope = FALSE)
    expect_length(imgs, 2L)
    expect_length(imgs[[1]], 2L)
    expect_equal(imgs[[1]][[1]], matrix(1:6, 2))
    expect_null(imgs[[1]][[2]])
    expect_null(client$get_ion_image_pixels("fixture", "C6H12O6", "+H", isotope_idx = 2))
})

test_that("METASPACE duplicate display masses and missing images preserve alignment", {
    client <- metaspace_review_client()
    local_mocked_bindings(metaspaceClient = function(...) client,
        .metaspace_image_pixels = function(img_info, ...) {
            if (is.null(img_info$url)) NULL else matrix(1:6, 2)
        })
    out <- getMetaspace("fixture", verbose = FALSE)
    expect_false(anyDuplicated(rownames(out$annotations)) > 0)
    expect_identical(names(out$images), rownames(out$annotations))
    expect_warning(mat <- metaspaceToFeatureMatrix(out, verbose = FALSE), "Omitting 1")
    expect_equal(nrow(mat), 1L)
    expect_identical(colnames(mat), c("1_1", "1_2", "2_1", "2_2", "3_1", "3_2"))
    expect_equal(as.numeric(mat), 1:6)
    expect_warning(spe <- loadMetaspace("fixture", verbose = FALSE), "Omitting 1")
    expect_equal(nrow(spe), 1L)
    expect_identical(SummarizedExperiment::rowData(spe)$mz, out$annotations$mz[1])
    expect_error(metaspaceToFeatureMatrix(list(images = list(NULL)), verbose = FALSE), "No downloadable")
})

test_that("METASPACE decoding supports grayscale, alpha masks and Python-style clipping", {
    pixels <- matrix(c(0, 0.2, 0.4, 0.6, 0.8, 1), 2)
    bytes <- png::writePNG(pixels, target = raw())
    response <- structure(list(status_code = 200L, headers = list(`content-type` = "image/png"),
        content = bytes, url = "https://example.invalid/image.png"), class = "response")
    local_mocked_bindings(GET = function(...) response, .package = "httr")
    metadata <- list(url = response$url, minIntensity = 0, maxIntensity = 100)
    image <- .metaspace_image_pixels(metadata, relative = FALSE)
    expect_equal(image, pixels * 100, tolerance = 1e-7)
    clipped <- .metaspace_image_pixels(metadata, relative = FALSE, hotspot_clipping = TRUE)
    expect_equal(max(clipped), unname(stats::quantile(c(20, 40, 60, 80, 100), 0.99)))
    rgba <- array(1, dim = c(1, 2, 4))
    rgba[1, 1, 4] <- 0
    response$content <- png::writePNG(rgba, target = raw())
    masked <- .metaspace_image_pixels(metadata, relative = FALSE)
    expect_identical(dim(masked), c(1L, 2L))
    expect_equal(as.numeric(masked), c(0, 100))
})

test_that("METASPACE authentication failure reports how to authenticate", {
    client <- metaspaceClient()
    environment(client$query)$POST <- function(...) structure(list(status_code = 401L), class = "response")
    expect_error(client$query("query { placeholder }"), "Authentication required")
})

test_that("rounded-mass collision suffixes never replace numeric feature masses", {
    client <- metaspace_review_client()
    local_mocked_bindings(metaspaceClient = function(...) client,
        .metaspace_image_pixels = function(...) matrix(1:6, 2))
    object <- loadMetaspace("fixture", verbose = FALSE)
    expect_identical(SummarizedExperiment::rowData(object)$mz, c(181.070664, 181.070665))
    expect_false(anyDuplicated(rownames(object)) > 0)
})
