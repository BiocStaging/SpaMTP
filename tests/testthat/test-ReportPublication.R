test_that("report publication supports an existing empty directory on Windows", {
    parent <- withr::local_tempdir()
    staging <- file.path(parent, "staging")
    destination <- file.path(parent, "report [1]")
    unrelated <- file.path(parent, "report 1")
    dir.create(staging)
    dir.create(destination)
    dir.create(unrelated)
    writeLines("complete report", file.path(staging, "report.html"))
    writeLines("keep", file.path(unrelated, "keep.txt"))
    publish <- getFromNamespace(".wf_publish_report", "SpaMTP")
    env <- new.env(parent = environment(publish))
    env$file.rename <- function(from, to) {
        # Model Windows: an existing destination directory blocks the rename.
        if (file.exists(to)) return(FALSE)
        base::file.rename(from, to)
    }
    environment(publish) <- env
    expect_identical(publish(staging, destination), destination)
    expect_identical(readLines(file.path(destination, "report.html")),
                     "complete report")
    expect_identical(readLines(file.path(unrelated, "keep.txt")), "keep")
    expect_false(dir.exists(staging))
})

test_that("files arriving during report rendering are preserved", {
    parent <- withr::local_tempdir()
    staging <- file.path(parent, "staging")
    destination <- file.path(parent, "report")
    dir.create(staging)
    dir.create(destination)
    writeLines("new report", file.path(staging, "report.html"))
    writeLines("user file", file.path(destination, ".keep"))
    publish <- getFromNamespace(".wf_publish_report", "SpaMTP")
    expect_error(publish(staging, destination), "existing files are preserved")
    expect_identical(readLines(file.path(destination, ".keep")), "user file")
    expect_true(file.exists(file.path(staging, "report.html")))
    unlink(destination, recursive = TRUE)
    writeLines("another file", destination)
    expect_error(publish(staging, destination), "existing files are preserved")
    expect_identical(readLines(destination), "another file")
})

test_that("report publication reports a failed rename and supports new paths", {
    parent <- withr::local_tempdir()
    staging <- file.path(parent, "staging")
    destination <- file.path(parent, "new report")
    dir.create(staging)
    writeLines("complete report", file.path(staging, "report.html"))
    publish <- getFromNamespace(".wf_publish_report", "SpaMTP")
    env <- new.env(parent = environment(publish))
    env$file.rename <- function(from, to) FALSE
    failedPublish <- publish
    environment(failedPublish) <- env
    expect_error(failedPublish(staging, destination), "Could not publish")
    expect_true(dir.exists(staging))
    publish(staging, destination)
    expect_true(file.exists(file.path(destination, "report.html")))
})
