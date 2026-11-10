# Development assistance: see inst/COPYRIGHTS.

# Event handlers record into a private environment. Combining rows once avoids
# repeated data-frame growth and does not modify the caller's environment.
.workflowRecords <- function(columns) {
    rows <- new.env(parent = emptyenv())
    list(
        add = function(...) {
            values <- list(...)
            if (length(values) != length(columns)) {
                stop(
                    "Incorrect number of workflow record fields.",
                    call. = FALSE
                )
            }
            key <- sprintf("%08d", length(rows) + 1L)
            assign(key, stats::setNames(values, columns), envir = rows)
            invisible(NULL)
        },
        count = function() length(rows),
        collect = function() {
            if (!length(rows)) {
                return(as.data.frame(
                    stats::setNames(
                        rep(list(character()), length(columns)),
                        columns
                    )
                ))
            }
            entries <- as.list(rows, sorted = TRUE)
            result <- do.call(rbind, lapply(entries, as.data.frame))
            rownames(result) <- NULL
            result
        }
    )
}
