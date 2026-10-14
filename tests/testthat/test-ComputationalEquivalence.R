test_that("blocked correlations agree with independent linear-model residuals", {
  withr::local_seed(872)
  n <- 40L
  group <- rep(c("A", "B"), each = n / 2L)
  trend <- seq_len(n)
  values <- matrix(rnorm(520L * n), 520L,
    dimnames = list(paste0("f", seq_len(520L)), paste0("p", seq_len(n))))
  values[1L, ] <- values[1L, ] + 2 * (group == "B")
  values[2L, ] <- 0
  values[3L, ] <- as.numeric(group == "B")
  values[4L, 1L] <- NA_real_
  values[257L, 2L] <- Inf
  values[513L, ] <- trend + 1e-8 * values[513L, ]
  design <- stats::model.matrix(~ group + trend)
  target <- values[1L, ]
  expected <- vapply(seq_len(nrow(values)), function(i) {
    y <- values[i, ]
    if (any(!is.finite(y)) || stats::sd(y) == 0) return(c(NA_real_, NA_real_))
    residuals <- stats::lm.fit(design, cbind(y, target))$residuals
    adjusted <- if (any(colSums(residuals^2) <
      1e-12 * c(sum((y - mean(y))^2), sum((target - mean(target))^2)))) {
      NA_real_
    } else stats::cor(residuals[, 1L], residuals[, 2L])
    c(adjusted, stats::cor(y, target))
  }, numeric(2))
  for (E in list(values, Matrix::Matrix(values, sparse = TRUE))) {
    object <- SpatialExperiment::SpatialExperiment(assays = list(counts = E),
      colData = S4Vectors::DataFrame(group = group, trend = trend,
        duplicate_trend = trend, row.names = colnames(E)),
      spatialCoords = cbind(x = trend, y = 0))
    result <- findCorrelatedFeatures(object, mz = "f1", nfeatures = NULL,
      covariates = c("group", "trend", "duplicate_trend"))
    ordered <- result[match(rownames(values), result$features), ]
    expect_equal(ordered$correlation_raw, expected[2L, ], tolerance = 1e-12)
    expect_equal(ordered$correlation, expected[1L, ], tolerance = 1e-12)
    expect_identical(attr(result, "association")$observations, colnames(values))
    expect_equal(attr(result, "association")$residual_df, n - qr(design)$rank)
    raw <- findCorrelatedFeatures(object, mz = "f1", nfeatures = NULL)
    expect_equal(raw$correlation, raw$correlation_raw)
    expect_equal(raw$correlation[match(rownames(values), raw$features)],
                 expected[2L, ], tolerance = 1e-12)
    for (reference in c("f2", "f3", "f4")) {
      invalid <- findCorrelatedFeatures(object, mz = reference, nfeatures = NULL,
        covariates = "group")
      expect_true(all(is.na(invalid$correlation)))
    }
    short <- findCorrelatedFeatures(object[, 1:2], mz = "f1", nfeatures = NULL)
    expect_true(all(is.na(short$correlation)))
  }
})

test_that("matrix-free Graph PCA retains the dense solution on wide inputs", {
  withr::local_seed(526)
  n <- 50L
  values <- matrix(rnorm(300L * n), 300L,
    dimnames = list(paste0("f", seq_len(300L)), paste0("p", seq_len(n))))
  for (E in list(values, Matrix::Matrix(values, sparse = TRUE))) {
    object <- SpatialExperiment::SpatialExperiment(assays = list(scaled = E),
      colData = S4Vectors::DataFrame(sample_id = rep(c("a", "b"), length.out = n),
        row.names = colnames(E)),
      spatialCoords = cbind(x = seq_len(n) %% 7, y = seq_len(n) %/% 7))
    fit <- runSpatialGraphPCA(object, n_components = 4, lambda = .7,
      platform = "ST", n_neighbors = 4, fast = TRUE, verbose = FALSE)
    A <- as.matrix(S4Vectors::metadata(fit)$spatialGraphs$SpatialKNN)
    input <- t(as.matrix(E))
    smoothed <- solve(diag(n) + .7 * (diag(rowSums(A)) - A), input)
    eig <- eigen(crossprod(input, smoothed), symmetric = TRUE)
    expected <- smoothed %*% eig$vectors[, 1:4]
    actual <- SingleCellExperiment::reducedDim(fit, "SpatialPCA")
    expect_equal(tcrossprod(actual), tcrossprod(expected),
      tolerance = 1e-7, ignore_attr = TRUE)
    expect_identical(SummarizedExperiment::assay(fit, "scaled"), E)
    expect_true(all(A[object$sample_id == "a", object$sample_id == "b"] == 0))
  }
})
