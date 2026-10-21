# Standalone installer for the Bioconductor submission. Sourcing this file
# only defines functions; call installSpaMTP() explicitly in a fresh R session.

installSpaMTP <- function(lib = NULL, buildVignettes = TRUE,
                         biocVersion = "devel", timeout = 600,
                         pandocDir = NULL, ref = "devel") {
    if (.spaMTPRVersion() < "4.6.0") {
        stop("SpaMTP requires R >= 4.6.0.", call. = FALSE)
    }
    if (!is.logical(buildVignettes) || length(buildVignettes) != 1L ||
        is.na(buildVignettes)) {
        stop("buildVignettes must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.numeric(timeout) || length(timeout) != 1L ||
        !is.finite(timeout) || timeout < 1) {
        stop("timeout must be a positive number of seconds.", call. = FALSE)
    }
    if (is.null(lib)) lib <- .spaMTPPersonalLibrary()
    for (value in list(lib, biocVersion, ref)) {
        if (!is.character(value) || length(value) != 1L ||
            is.na(value) || !nzchar(value)) {
            stop("lib, biocVersion and ref must be nonempty strings.",
                 call. = FALSE)
        }
    }
    if (!is.null(pandocDir) && (!is.character(pandocDir) ||
        length(pandocDir) != 1L || is.na(pandocDir) ||
        !dir.exists(pandocDir))) {
        stop("pandocDir must name an existing directory containing Pandoc.",
             call. = FALSE)
    }
    if ("SpaMTP" %in% loadedNamespaces()) {
        stop("Restart R before installing SpaMTP, especially on Windows.",
             call. = FALSE)
    }

    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)
    if (file.access(lib, 2L) != 0L) {
        stop("Choose a writable personal package library with lib=.",
             call. = FALSE)
    }
    previousLibs <- .libPaths()
    previousOptions <- options("timeout", "repos")
    envNames <- c("R_LIBS", "R_LIBS_USER", "R_DEFAULT_INTERNET_TIMEOUT",
                  "RSTUDIO_PANDOC")
    previousEnv <- Sys.getenv(envNames, unset = NA_character_)
    on.exit({
        .libPaths(previousLibs)
        options(previousOptions)
        .spaMTPRestoreInstallEnv(previousEnv)
    }, add = TRUE)
    .libPaths(c(lib, previousLibs))
    inheritedTimeout <- suppressWarnings(as.numeric(c(
        getOption("timeout"), Sys.getenv("R_DEFAULT_INTERNET_TIMEOUT"))))
    timeout <- max(timeout, inheritedTimeout[is.finite(inheritedTimeout)])
    options(timeout = timeout)
    libraryPaths <- paste(.libPaths(), collapse = .Platform$path.sep)
    Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = as.character(ceiling(timeout)),
               R_LIBS = libraryPaths, R_LIBS_USER = libraryPaths)
    repos <- getOption("repos")
    if (is.null(names(repos)) || !"CRAN" %in% names(repos) ||
        is.na(repos[["CRAN"]]) || !nzchar(repos[["CRAN"]]) ||
        repos[["CRAN"]] == "@CRAN@") {
        repos <- c(repos[names(repos) != "CRAN"],
                   CRAN = "https://cloud.r-project.org")
    }
    options(repos = repos)
    message("Download timeout: ", timeout, " seconds; library: ", lib)

    ops <- .spaMTPInstallOperations()
    bootstrap <- c("BiocManager", "remotes",
                   if (buildVignettes) "rmarkdown")
    for (package in bootstrap) {
        if (!ops$installed(package, "0")) {
            ops$bootstrap(package, lib, repos[["CRAN"]])
        }
        ops$check(package, "0")
    }
    # Fail before the large dependency installation if R cannot build articles.
    if (buildVignettes) {
        pandoc <- ops$pandoc(pandocDir)
        Sys.setenv(RSTUDIO_PANDOC = pandoc$dir)
        message("Pandoc ", pandoc$version, ": ", pandoc$dir)
    }
    ops$bioc(version = biocVersion, lib = lib, ask = FALSE, update = FALSE)

    # Data has graduated to Bioconductor. The database may still be in review;
    # prefer an adequate repository version whenever it becomes available.
    if (!ops$installed("SpaMTPdb", "0.99.4")) {
        available <- ops$available()
        official <- "SpaMTPdb" %in% rownames(available) &&
            numeric_version(available["SpaMTPdb", "Version"]) >= "0.99.4"
        target <- if (official) "SpaMTPdb" else "BiocStaging/SpaMTPdb@devel"
        message("Installing database dependency from ", target)
        ops$bioc(target, lib = lib, ask = FALSE, update = FALSE,
                 build_vignettes = FALSE)
    }
    ops$check("SpaMTPdb", "0.99.4")
    if (!ops$installed("SpaMTPData", "0.99.5")) {
        ops$bioc("SpaMTPData", lib = lib, ask = FALSE, update = FALSE)
    }
    ops$check("SpaMTPData", "0.99.5")
    ops$check("magick", "0")
    ops$check("SpatialExperiment", "0")
    ops$bioc(paste0("BiocStaging/SpaMTP@", ref), lib = lib,
             dependencies = TRUE, build_vignettes = buildVignettes,
             force = buildVignettes, ask = FALSE, update = FALSE)
    ops$check("SpaMTP", "0.99.11")
    message("SpaMTP installation and namespace checks completed.")
    invisible(lib)
}

# Separate from the guard so tests can exercise the installer on any R.
.spaMTPRVersion <- function() getRversion()

.spaMTPPersonalLibrary <- function() {
    personal <- Sys.getenv("R_LIBS_USER")
    if (nzchar(personal)) {
        paths <- strsplit(personal, .Platform$path.sep, fixed = TRUE)[[1L]]
        paths <- paths[nzchar(paths)]
        if (length(paths)) return(path.expand(paths[[1L]]))
    }
    stop("Set lib= to a writable personal package directory, or set R_LIBS_USER.",
         call. = FALSE)
}

.spaMTPRestoreInstallEnv <- function(values) {
    absent <- is.na(values)
    if (any(absent)) Sys.unsetenv(names(values)[absent])
    if (any(!absent)) do.call(Sys.setenv, as.list(values[!absent]))
    invisible(NULL)
}

.spaMTPInstallOperations <- function() {
    list(
        installed = function(package, version) {
            tryCatch(utils::packageVersion(package) >= version,
                     error = function(e) FALSE)
        },
        bootstrap = function(package, lib, repos) {
            utils::install.packages(package, lib = lib, repos = repos)
        },
        bioc = function(...) BiocManager::install(...),
        available = function() {
            utils::available.packages(repos = BiocManager::repositories(),
                                      type = "source")
        },
        check = function(package, version) {
            tryCatch({
                ns <- loadNamespace(package)
                if (getNamespaceVersion(ns) < version) {
                    stop("installed version is older than ", version)
                }
            }, error = function(e) {
                stop("Cannot load ", package, ": ", conditionMessage(e),
                     ". Inspect the preceding download/install messages; ",
                     "a 60-second ZIP timeout requires a larger timeout. ",
                     "Restart R after repairing an incomplete installation.",
                     call. = FALSE)
            })
            invisible(TRUE)
        },
        pandoc = .spaMTPInstallPandoc
    )
}

.spaMTPWindowsPandocDirs <- function(
        localAppData = Sys.getenv("LOCALAPPDATA"),
        programFiles = Sys.getenv("ProgramFiles"),
        programFilesX86 = Sys.getenv("ProgramFiles(x86)")) {
    roots <- c(if (nzchar(localAppData)) file.path(localAppData, "Programs"),
               programFiles[nzchar(programFiles)],
               programFilesX86[nzchar(programFilesX86)])
    layouts <- c("resources/app/bin/quarto/bin/tools", "bin/quarto/bin/tools",
                 "resources/app/bin/pandoc", "bin/pandoc")
    unique(c(unlist(lapply(roots, function(root) {
        file.path(root, "RStudio", layouts)
    }), use.names = FALSE),
    if (nzchar(localAppData)) file.path(localAppData, "Pandoc"),
    file.path(roots, "Pandoc")))
}

.spaMTPInstallPandoc <- function(pandocDir = NULL) {
    if (!is.null(pandocDir)) {
        found <- rmarkdown::find_pandoc(cache = FALSE, dir = pandocDir)
    } else {
        found <- rmarkdown::find_pandoc(cache = FALSE)
        if (found$version < "1.12.3" && .Platform$OS.type == "windows") {
            candidates <- .spaMTPWindowsPandocDirs()
            candidates <- candidates[file.exists(file.path(candidates, "pandoc.exe"))]
            if (length(candidates)) {
                found <- rmarkdown::find_pandoc(cache = FALSE, dir = candidates)
            }
        }
    }
    if (is.null(found$dir) || found$version < "1.12.3") {
        stop("Pandoc is required to build vignettes. Run from RStudio, ",
             "install Pandoc, or pass pandocDir= with its containing directory. ",
             "Set buildVignettes=FALSE only if you do not need local articles.",
             call. = FALSE)
    }
    found$dir <- normalizePath(found$dir, winslash = "/", mustWork = TRUE)
    found
}
