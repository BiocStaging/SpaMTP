installerEnvironment <- function() {
    path <- testthat::test_path("..", "..", "inst", "workflows", "install_spamtp.R")
    if (!file.exists(path)) {
        path <- system.file("workflows", "install_spamtp.R", package = "SpaMTP")
    }
    env <- new.env(parent = globalenv())
    sys.source(path, envir = env)
    # Installer logic is tested independently of the R running the tests.
    env$.spaMTPRVersion <- function() numeric_version("4.6.0")
    env
}

fakeInstaller <- function(dbVersion = NULL, failCheck = NULL) {
    env <- installerEnvironment()
    seen <- new.env(parent = emptyenv())
    seen$calls <- list()
    seen$checks <- character()
    record <- function(...) {
        seen$calls[[length(seen$calls) + 1L]] <- list(...)
    }
    # A new R process has not loaded SpaMTP; the test runner necessarily has.
    env$loadedNamespaces <- function() character()
    env$.spaMTPInstallOperations <- function() list(
        installed = function(package, version) FALSE,
        bootstrap = function(package, lib, repos) {
            record(kind = "bootstrap", package = package, lib = lib, repos = repos)
        },
        bioc = function(...) {
            record(kind = "bioc", args = list(...), timeout = getOption("timeout"),
                   childTimeout = Sys.getenv("R_DEFAULT_INTERNET_TIMEOUT"),
                   childLibs = Sys.getenv("R_LIBS_USER"),
                   inheritedLibs = Sys.getenv("R_LIBS"),
                   pandoc = Sys.getenv("RSTUDIO_PANDOC"))
        },
        available = function() {
            if (is.null(dbVersion)) {
                matrix(character(), 0, 1, dimnames = list(NULL, "Version"))
            } else {
                matrix(dbVersion, 1, 1, dimnames = list("SpaMTPdb", "Version"))
            }
        },
        check = function(package, version) {
            seen$checks <- c(seen$checks, package)
            if (identical(package, failCheck)) stop("Dependency load failed")
        },
        pandoc = function(dir) {
            record(kind = "pandoc")
            list(dir = dir, version = numeric_version("3.10"))
        }
    )
    list(env = env, seen = seen)
}

test_that("the installer can be sourced without installing or loading anything", {
    before <- loadedNamespaces()
    env <- installerEnvironment()
    expect_setequal(loadedNamespaces(), before)
    expect_true(all(vapply(as.list(env), is.function, logical(1))))
    expect_error(env$installSpaMTP(timeout = 0), "positive")
    expect_error(env$installSpaMTP(buildVignettes = NA), "TRUE or FALSE")
    if ("SpaMTP" %in% loadedNamespaces()) {
        expect_error(env$installSpaMTP(), "Restart R")
    }
})

test_that("the installer stops on R older than 4.6.0", {
    env <- installerEnvironment()
    env$.spaMTPRVersion <- function() numeric_version("4.5.3")
    expect_error(env$installSpaMTP(), "requires R >= 4.6.0")
})

test_that("installation handles an empty library and a database still in staging", {
    fake <- fakeInstaller()
    lib <- withr::local_tempdir(pattern = "personal library ")
    pandoc <- withr::local_tempdir(pattern = "Pandoc tools ")
    withr::local_options(list(timeout = 60, repos = c(CRAN = "@CRAN@")))
    withr::local_envvar(c(R_DEFAULT_INTERNET_TIMEOUT = NA, RSTUDIO_PANDOC = NA))
    before <- list(options = options("timeout", "repos"), libs = .libPaths(),
        env = Sys.getenv(c("R_LIBS", "R_LIBS_USER", "R_DEFAULT_INTERNET_TIMEOUT",
                           "RSTUDIO_PANDOC"), unset = NA_character_))
    fake$env$installSpaMTP(lib = lib, pandocDir = pandoc, biocVersion = "3.24")
    calls <- fake$seen$calls
    bootstrap <- Filter(function(x) x$kind == "bootstrap", calls)
    expect_identical(vapply(bootstrap, `[[`, "", "package"),
                     c("BiocManager", "remotes", "rmarkdown"))
    bioc <- Filter(function(x) x$kind == "bioc", calls)
    expect_identical(bioc[[1]]$args$version, "3.24")
    expect_identical(bioc[[2]]$args[[1]], "BiocStaging/SpaMTPdb@devel")
    expect_identical(bioc[[3]]$args[[1]], "SpaMTPData")
    expect_identical(bioc[[4]]$args[[1]], "BiocStaging/SpaMTP@devel")
    expect_true(bioc[[4]]$args$build_vignettes)
    expect_true(bioc[[4]]$args$force)
    expect_true(bioc[[4]]$args$dependencies)
    expect_true(all(vapply(bioc, function(x) x$timeout == 600 &&
        x$childTimeout == "600" && x$pandoc == pandoc, logical(1))))
    expect_identical(strsplit(bioc[[4]]$childLibs, .Platform$path.sep,
                             fixed = TRUE)[[1]][1], normalizePath(lib, winslash = "/"))
    expect_identical(bioc[[4]]$inheritedLibs, bioc[[4]]$childLibs)
    expect_identical(options("timeout", "repos"), before$options)
    expect_identical(.libPaths(), before$libs)
    expect_identical(Sys.getenv(names(before$env), unset = NA_character_), before$env)
})

test_that("official database versions are preferred only when sufficient", {
    for (version in c("0.99.3", "0.99.4", "1.0.0")) {
        fake <- fakeInstaller(dbVersion = version)
        fake$env$installSpaMTP(lib = withr::local_tempdir(), buildVignettes = FALSE)
        calls <- Filter(function(x) x$kind == "bioc", fake$seen$calls)
        expected <- if (version == "0.99.3") "BiocStaging/SpaMTPdb@devel" else "SpaMTPdb"
        expect_identical(calls[[2]]$args[[1]], expected)
        expect_false(calls[[4]]$args$build_vignettes)
        expect_false(calls[[4]]$args$force)
        expect_false(any(vapply(fake$seen$calls, function(x) x$kind == "pandoc",
                               logical(1))))
    }
})

test_that("installation preserves larger timeouts and cleans up on dependency failure", {
    fake <- fakeInstaller(failCheck = "SpaMTPData")
    withr::local_options(timeout = 900)
    withr::local_envvar(c(R_LIBS = "previous library",
                         R_DEFAULT_INTERNET_TIMEOUT = "1200",
                         RSTUDIO_PANDOC = "previous tools"))
    before <- .libPaths()
    expect_error(fake$env$installSpaMTP(lib = withr::local_tempdir(),
                 buildVignettes = FALSE), "Dependency load failed")
    calls <- Filter(function(x) x$kind == "bioc", fake$seen$calls)
    expect_length(calls, 3L)
    expect_true(all(vapply(calls, function(x) x$timeout == 1200 &&
        x$childTimeout == "1200", logical(1))))
    expect_equal(getOption("timeout"), 900)
    expect_identical(Sys.getenv("R_DEFAULT_INTERNET_TIMEOUT"), "1200")
    expect_identical(Sys.getenv("R_LIBS"), "previous library")
    expect_identical(Sys.getenv("RSTUDIO_PANDOC"), "previous tools")
    expect_identical(.libPaths(), before)
})

test_that("missing Pandoc stops before the full dependency installation", {
    fake <- fakeInstaller()
    ops <- fake$env$.spaMTPInstallOperations()
    ops$pandoc <- function(dir) stop("Pandoc unavailable")
    fake$env$.spaMTPInstallOperations <- function() ops
    withr::local_options(timeout = 60)
    expect_error(fake$env$installSpaMTP(lib = withr::local_tempdir()),
                 "Pandoc unavailable")
    expect_false(any(vapply(fake$seen$calls, function(x) x$kind == "bioc",
                           logical(1))))
    expect_equal(getOption("timeout"), 60)
})

test_that("Windows Pandoc discovery covers user and system RStudio layouts", {
    env <- installerEnvironment()
    dirs <- env$.spaMTPWindowsPandocDirs("C:/Users/Test/AppData/Local",
                                       "C:/Program Files", "")
    expect_true("C:/Users/Test/AppData/Local/Programs/RStudio/resources/app/bin/quarto/bin/tools" %in% dirs)
    expect_true("C:/Program Files/RStudio/bin/quarto/bin/tools" %in% dirs)
    expect_true("C:/Program Files/Pandoc" %in% dirs)
    expect_true("C:/Users/Test/AppData/Local/Pandoc" %in% dirs)
    expect_length(env$.spaMTPWindowsPandocDirs("", "", ""), 0L)
})

test_that("the default personal library need not already exist", {
    env <- installerEnvironment()
    personal <- file.path(withr::local_tempdir(), "new personal library")
    withr::local_envvar(c(R_LIBS_USER = personal))
    expect_identical(env$.spaMTPPersonalLibrary(), personal)
    fake <- fakeInstaller()
    fake$env$installSpaMTP(buildVignettes = FALSE)
    expect_true(dir.exists(personal))
    calls <- Filter(function(x) x$kind == "bioc", fake$seen$calls)
    expect_identical(calls[[1]]$args$lib, normalizePath(personal, winslash = "/"))
    Sys.unsetenv("R_LIBS_USER")
    expect_error(env$.spaMTPPersonalLibrary(), "writable personal")
})

test_that("fresh build processes inherit the timeout, library and Pandoc directory", {
    fake <- fakeInstaller()
    lib <- withr::local_tempdir(pattern = "child library ")
    inherited <- withr::local_tempdir(pattern = "inherited library ")
    pandoc <- withr::local_tempdir(pattern = "child Pandoc ")
    output <- tempfile(fileext = ".rds")
    script <- tempfile(fileext = ".R")
    withr::defer(unlink(c(script, output)))
    writeLines(c(
        "saveRDS(list(timeout = getOption('timeout'), libs = .libPaths(),",
        "pandoc = Sys.getenv('RSTUDIO_PANDOC')), commandArgs(TRUE)[1])"
    ), script)
    ops <- fake$env$.spaMTPInstallOperations()
    ops$bioc <- function(...) {
        if (is.null(list(...)$version)) return(invisible(NULL))
        executable <- file.path(R.home("bin"),
            if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
        status <- system2(executable,
            c("--vanilla", shQuote(script), shQuote(output)),
            stdout = TRUE, stderr = TRUE)
        expect_null(attr(status, "status"))
    }
    fake$env$.spaMTPInstallOperations <- function() ops
    withr::local_envvar(c(R_LIBS = inherited, R_DEFAULT_INTERNET_TIMEOUT = NA))
    withr::local_options(timeout = 60)
    fake$env$installSpaMTP(lib = lib, timeout = 725, pandocDir = pandoc)
    child <- readRDS(output)
    expect_equal(child$timeout, 725)
    expect_identical(child$libs[[1L]], normalizePath(lib, winslash = "/"))
    expect_identical(child$pandoc, pandoc)
    expect_identical(Sys.getenv("R_LIBS"), inherited)
})

test_that("Pandoc discovery validates executables, not just directories", {
    skip_if_not_installed("rmarkdown")
    env <- installerEnvironment()
    empty <- withr::local_tempdir()
    # Recompute discovery after the test so an absent explicit path is not cached.
    withr::defer(rmarkdown::find_pandoc(cache = FALSE))
    expect_error(env$.spaMTPInstallPandoc(empty), "Pandoc is required")
    normal <- rmarkdown::find_pandoc(cache = FALSE)
    if (!is.null(normal$dir) && normal$version >= "1.12.3") {
        found <- env$.spaMTPInstallPandoc(normal$dir)
        expect_identical(found$version, normal$version)
        expect_identical(found$dir, normalizePath(normal$dir, winslash = "/"))
    }
})
