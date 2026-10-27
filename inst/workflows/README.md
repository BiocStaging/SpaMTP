# SpaMTP installation and workflow recipes

SpaMTP and SpaMTPdb are under Bioconductor review. SpaMTPData is available
from Bioconductor. Install the data packages before the analysis package.
The packages require R >= 4.6.0; the commands here select Bioconductor 3.24.

## Standalone installer

Start a fresh R session, particularly on Windows, where loaded packages can
prevent replacement of installed files.

```r
source("https://raw.githubusercontent.com/BiocStaging/SpaMTP/devel/inst/workflows/install_spamtp.R")
installSpaMTP(biocVersion = "3.24")
```

Sourcing `install_spamtp.R` only defines functions. `installSpaMTP()` checks R,
raises the download timeout to at least 600 seconds, installs BiocManager and
remotes if needed, and checks Pandoc before installing the larger dependency
set. It installs SpaMTPdb (>= 0.99.4), SpaMTPData (>= 0.99.5), then SpaMTP.
It prefers an adequate SpaMTPdb version from the configured Bioconductor
repositories and otherwise uses BiocStaging. SpaMTP itself comes from
`BiocStaging/SpaMTP`, branch `devel`.

The installer checks that the packages load, restores session options and
environment variables on exit, and gives child R processes the selected
library first. It preserves an existing longer timeout.

| Argument | Default | Use |
|---|---|---|
| `lib` | First directory in `R_LIBS_USER`, created if needed | Choose a writable library; add it to `.libPaths()` in later sessions. |
| `buildVignettes` | `TRUE` | Build all seven articles. `FALSE` skips the articles and Pandoc check. |
| `biocVersion` | `"devel"` | Select a Bioconductor version compatible with your R installation; the example pins `"3.24"`. |
| `timeout` | `600` | Minimum download timeout in seconds. |
| `ref` | `"devel"` | Select a branch or fixed commit SHA of BiocStaging/SpaMTP. |
| `pandocDir` | `NULL` | Directory containing `pandoc` or `pandoc.exe`. |

With `buildVignettes = TRUE`, SpaMTP is rebuilt even if the same commit was
previously installed without articles. For slower connections or a separate
Pandoc installation:

```r
installSpaMTP(biocVersion = "3.24", timeout = 1800,
              pandocDir = "C:/Program Files/Pandoc")
```

## Manual installation during review

Run in RStudio, or make Pandoc visible to R before building the vignettes.

```r
options(timeout = max(600, getOption("timeout")))
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = as.character(getOption("timeout")))
install.packages(c("BiocManager", "remotes", "rmarkdown"),
                 repos = "https://cloud.r-project.org")
BiocManager::install(version = "3.24", ask = FALSE, update = FALSE)
stopifnot(rmarkdown::pandoc_available())

BiocManager::install("BiocStaging/SpaMTPdb@devel", ask = FALSE, update = FALSE)
BiocManager::install("SpaMTPData", ask = FALSE, update = FALSE)
loadNamespace("SpaMTPData")
BiocManager::install("BiocStaging/SpaMTP@devel", dependencies = TRUE,
                     build_vignettes = TRUE, ask = FALSE, update = FALSE)
```

These manual commands change the current session settings. After both staging
packages are accepted and available in your Bioconductor release,
`BiocManager::install("SpaMTP")` will install the analysis package and its
required dependencies. Install SpaMTPData separately for the example workflows.

## Local source checkouts

Install the dependencies listed in each `DESCRIPTION` first. From the SpaMTP
checkout, with the companion checkouts alongside it:

```sh
R CMD INSTALL ../SpaMTPdb
R CMD INSTALL ../SpaMTPData
R CMD INSTALL .
```

`R CMD INSTALL .` installs the source but does not build missing vignettes.
For an archive containing the articles, run `R CMD build .` and install the
resulting tarball. Keep generated archives and check output outside the source
checkout, or in its ignored `.local/` directory.

## Troubleshooting

* **Missing SpaMTPdb:** while it is under review, install it from BiocStaging
  before SpaMTP. SpaMTPData is installed from Bioconductor.
* **Missing remotes:** BiocManager uses remotes to install GitHub packages;
  install it before requesting a staging package.
* **Missing Pandoc:** run in RStudio, pass `pandocDir`, or set `RSTUDIO_PANDOC`
  to the directory containing the executable. The installer searches common
  RStudio locations on Windows. `buildVignettes = FALSE` skips article builds.
* **Missing magick on Windows:** inspect earlier download messages. The
  reported fresh-library failure followed a dependency download timeout;
  increase `timeout` and retry. This does not establish the cause of every
  magick installation failure. The Windows binary normally includes the
  required libraries; see the [magick installation guide](https://docs.ropensci.org/magick/#installation).
* **Repository-index HTTP 404:** check whether the failed index is for Windows
  binaries of a data package that is supplied as source. Inspect the actual
  package-installation result before treating the index warning as a failure.
  Retain binary support for compiled dependencies, and use matching Rtools
  when source compilation is required.

## Installed workflow recipes

| File | Purpose |
|---|---|
| `install_spamtp.R` | Install the submission and its dependencies. |
| `prepare_mouse_brain.R` | Prepare native mouse-brain inputs with explicit spatial configuration. |
| `run_mouse_brain.R` | Run the analysis and produce an HTML report with exported tables. |
| `prepare_vignette_data.R` | Reproduce the compact measured datasets used by the introduction. |

Inspect each script's usage before running it. The numbered SpaMTP vignettes
explain the analysis choices and resource versions. Find these recipes after
installation with `system.file("workflows", package = "SpaMTP")`.
