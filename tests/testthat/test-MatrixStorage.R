# Development assistance: see inst/COPYRIGHTS.

test_that("CSV import chooses storage without changing values or metadata", {
    cases <- list(
        dense = rbind(c(1, 2, 3), c(4, 5, 6)),
        sparse = rbind(c(0, 2, 0), c(0, 0, 6)),
        half = rbind(c(1, 0, 3), c(0, 5, 0)),
        zero = matrix(0, 2, 3),
        missing = rbind(c(NA, 0, 0), c(0, 0, 6))
    )
    for (name in names(cases)) {
        values <- cases[[name]]
        table <- data.frame(
            x = c(1, 2, 1), y = c(1, 1, 2),
            mz100.1 = values[1, ], mz200.2 = values[2, ]
        )
        path <- withr::local_tempfile(fileext = ".csv")
        utils::write.csv(table, path, row.names = FALSE)
        object <- readSMMatrix(path,
            mz.prefix = "mz", project.name = "sample", verbose = FALSE
        )
        counts <- SummarizedExperiment::assay(object, "counts")
        if (name == "dense") {
            expect_true(is.matrix(counts))
            expect_type(counts, "double")
        } else {
            expect_s4_class(counts, "dgCMatrix")
        }
        expect_equal(unname(as.matrix(counts)), values)
        expect_identical(rownames(object), c("mz-100.1", "mz-200.2"))
        expect_identical(colnames(object), c("1_1", "2_1", "1_2"))
        expect_equal(SummarizedExperiment::rowData(object)$mz, c(100.1, 200.2))
        expect_equal(
            unname(SpatialExperiment::spatialCoords(object)),
            unname(as.matrix(table[c("x", "y")]))
        )
        expect_identical(object$orig.ident, rep("sample", 3))
        expect_identical(SingleCellExperiment::mainExpName(object), "Spatial")
        expect_identical(
            S4Vectors::metadata(object)$spamtp_assay_name, "Spatial"
        )
    }
})

test_that("single-feature CSVs retain a two-dimensional assay", {
    path <- withr::local_tempfile(fileext = ".csv")
    table <- data.frame(
        ignored = letters[1:3], x = 1:3, y = 1,
        mz100.1 = c(0, 7, 0)
    )
    utils::write.csv(table, path, row.names = FALSE)
    object <- readSMMatrix(path,
        feature.start.column = 2, mz.prefix = "mz", verbose = FALSE
    )
    expect_identical(dim(object), c(1L, 3L))
    counts <- SummarizedExperiment::assay(object, "counts")
    expect_s4_class(counts, "dgCMatrix")
    expect_equal(unname(as.matrix(counts)), matrix(c(0, 7, 0), nrow = 1))
    expect_identical(rownames(object), "mz-100.1")
    expect_identical(colnames(object), c("1_1", "2_1", "3_1"))
})

test_that("square sparse CSVs retain all values in a general sparse matrix", {
    path <- withr::local_tempfile(fileext = ".csv")
    table <- data.frame(
        x = 1:3, y = 1,
        mz100 = c(1, 0, 0), mz200 = c(2, 3, 0), mz300 = c(0, 0, 4)
    )
    utils::write.csv(table, path, row.names = FALSE)
    object <- readSMMatrix(path, mz.prefix = "mz", verbose = FALSE)
    counts <- SummarizedExperiment::assay(object, "counts")
    expect_s4_class(counts, "dgCMatrix")
    expect_equal(unname(as.matrix(counts)), unname(t(as.matrix(table[3:5]))))
})

test_that("large integer intensities are read as double values", {
    path <- withr::local_tempfile(fileext = ".csv")
    writeLines(c("x,y,mz100", "1,1,3000000000", "2,1,4000000000"), path)
    object <- readSMMatrix(path, mz.prefix = "mz", verbose = FALSE)
    counts <- SummarizedExperiment::assay(object, "counts")
    expect_type(counts, "double")
    expect_equal(unname(counts), matrix(c(3e9, 4e9), nrow = 1))
})

test_that("CSV import diagnoses missing and nonnumeric intensity columns", {
    path <- withr::local_tempfile(fileext = ".csv")
    writeLines(c("x,y", "1,1", "2,1"), path)
    expect_error(readSMMatrix(path, verbose = FALSE), "at least one m/z")
    writeLines(c("x,y,mz100", "1,1,text", "2,1,3"), path)
    expect_error(
        readSMMatrix(path, mz.prefix = "mz", verbose = FALSE),
        "Intensity columns must contain numeric"
    )
})

test_that("compaction preserves delayed assays and default dense storage", {
    values <- matrix(c(0, 0, 0, 1), 2)
    expect_identical(.compactAssay(values), values)
    delayed <- DelayedArray::DelayedArray(values)
    expect_identical(.compactAssay(delayed, sparsify = TRUE), delayed)
    denseMissing <- matrix(c(NA, NA, 1, 0), 2)
    expect_identical(.compactAssay(denseMissing, sparsify = TRUE), denseMissing)
})

test_that("normalization rejects invalid sparse and delayed intensities", {
    for (value in c(Inf, -Inf, NA_real_, NaN, -1)) {
        sparse <- Matrix::sparseMatrix(
            i = c(1, 2), j = c(1, 4), x = c(value, 1), dims = c(2, 4)
        )
        for (counts in list(sparse, DelayedArray::DelayedArray(sparse))) {
            object <- SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = counts),
                colData = S4Vectors::DataFrame(group = c("A", "A", "B", "B"))
            )
            expect_error(normalizeSMData(object, verbose = FALSE),
                "finite, non-negative"
            )
            expect_error(tmmNormalize(object, ident = "group", refIdent = "A",
                verbose = FALSE
            ), "finite, non-negative")
        }
    }
})
