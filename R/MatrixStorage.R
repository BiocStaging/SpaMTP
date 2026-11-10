# Development assistance: see inst/COPYRIGHTS.

# Matrix storage helpers for import and preprocessing.
#
# MSI intensity matrices are often nearly dense, so a dgCMatrix can cost more
# memory than a plain matrix: it stores an 8-byte value plus a 4-byte row index
# for every nonzero entry. Results are therefore kept sparse only while at most
# half of their entries are nonzero. Import can also compact a dense input;
# preprocessing preserves its input representation unless sparse storage is
# wasteful. DelayedArray input stays delayed.

.isDelayed <- function(x) methods::is(x, "DelayedArray")

.compactAssay <- function(x, maxDensity = 0.5, sparsify = FALSE) {
    if (.isDelayed(x) || !length(x)) {
        return(x)
    }
    sparse <- methods::is(x, "sparseMatrix")
    if (!sparse && !sparsify) {
        return(x)
    }
    # Missing entries also occupy space in a sparse matrix.
    density <- Matrix::nnzero(x, na.counted = TRUE) /
        prod(as.numeric(dim(x)))
    if (density > maxDensity) {
        return(if (sparse) as.matrix(x) else x)
    }
    if (sparse) x else methods::as(
        methods::as(x, "CsparseMatrix"), "generalMatrix"
    )
}

# Divide each column by one factor without densifying sparse input.
.divideColumns <- function(x, factors) {
    if (.isDelayed(x)) {
        return(DelayedArray::sweep(x, 2L, factors, "/"))
    }
    if (methods::is(x, "sparseMatrix")) {
        return(x %*% Matrix::Diagonal(x = 1 / factors))
    }
    x / rep(factors, each = nrow(x))
}

.columnTotals <- function(x) {
    if (.isDelayed(x)) DelayedArray::colSums(x) else Matrix::colSums(x)
}

# Testing !is.finite(sparse) materializes implicit zeros as TRUE. Testing NA
# and infinity separately keeps those zeros implicit (and delayed inputs lazy).
.allFiniteAssay <- function(x) !anyNA(x) && !any(is.infinite(x))
