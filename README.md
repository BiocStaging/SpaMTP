# SpaMTP <img src="man/figures/logo.png" align="right" height="100" alt="SpaMTP logo" />

**SpaMTP** is an R package for analysing spatial metabolomics and integrating
it with spatial transcriptomics. It connects mass spectrometry imaging (MSI)
preprocessing, candidate metabolite annotation, spatial statistics, pathway
analysis and interactive reports.

The workflow uses Bioconductor containers: Cardinal for raw MSI spectra and
`SpatialExperiment` for aligned or binned measurements, spatial coordinates
and paired modalities. Seurat conversion is optional and explicit.

## Installation

SpaMTP and SpaMTPdb are under Bioconductor review. The source packages require
**R >= 4.6.0**. The installation example uses **Bioconductor 3.24**; Pandoc
(included with RStudio) is needed to build the vignettes.

Start a fresh R session and run:

```r
source("https://raw.githubusercontent.com/BiocStaging/SpaMTP/devel/inst/workflows/install_spamtp.R")
installSpaMTP(biocVersion = "3.24")
```

The installer raises the download timeout, installs the companion packages
before SpaMTP, and builds all seven vignettes. It uses a personal R library
and restores session settings afterwards. For an installation without local
vignettes, use `buildVignettes = FALSE`.

| Package | Role | Installation source during review |
|---|---|---|
| [SpaMTPdb](https://github.com/BCRL-tylu/SpaMTPdb) >= 0.99.4 | Versioned RaMP annotation resources and a human HGNC reference | BiocStaging `devel`, or Bioconductor when available |
| [SpaMTPData](https://bioconductor.org/packages/SpaMTPData/) >= 0.99.5 | Experimental datasets and example data | Bioconductor |
| SpaMTP | Analysis, annotation, visualisation and integration | BiocStaging `devel` |

See the [installation and workflow guide](inst/workflows/README.md) for
manual installation, local checkouts, Windows/Pandoc setup and installer
options. SpaMTPdb must be installed explicitly before a manual staging
installation of SpaMTP.

## Quick start

This example uses a compact measured mouse-bladder MSI panel bundled with
SpaMTP, so it needs no data download after installation. The 79 features and
2,210 retained locations preserve the source panel and a regular tissue grid.

```r
library(SpaMTP)

bladder <- readRDS(system.file(
    "extdata", "bladder_demo.rds", package = "SpaMTP"
))
dim(bladder)
SummarizedExperiment::assayNames(bladder)

# Display one measured ion using its archived, processed intensities.
plotSpatialFeature(bladder, features = rownames(bladder)[1L],
                   assayName = "counts")
```

The [introductory vignette](vignettes/SpaMTP_Bioconductor.Rmd) explains the
tissue compartments, ion selection and processing history. These are
single-section examples: spatial differences describe the specimen, and
accurate-mass matches are candidate identities that need further chemical
evidence.

## Tutorials

Read the vignettes in the numbered order, or start with the topic needed for
your analysis. After installing with vignettes, use `browseVignettes("SpaMTP")`.

| Order | Vignette | Covers |
|---|---|---|
| 1 | [Spatial metabolomics with native containers](vignettes/SpaMTP_Bioconductor.Rmd) | Small mouse-bladder and pig-tissue examples, with biological context |
| 2 | [Data access and migration](vignettes/Modern_Data_Access.Rmd) | Imports, metadata, coordinates and optional Seurat conversion |
| 3 | [Companion data and annotation resources](vignettes/SpaMTP_Companion_Workflow.Rmd) | SpaMTPData experiments and verified SpaMTPdb resources |
| 4 | [Metabolite annotation and pathway follow-up](vignettes/Metabolite_Annotation_Pipeline.Rmd) | Accurate-mass candidates, structural evidence and pathway interpretation |
| 5 | [Gene mapping and joint pathway analysis](vignettes/Gene_Identifier_Mapping.Rmd) | Human gene identifiers, measured backgrounds and shared pathways |
| 6 | [Pathway resources and network provenance](vignettes/Pathway_Database_Integration.Rmd) | Membership, topology and resource provenance |
| 7 | [Spatial multi-omics analysis and evaluation](vignettes/End_to_End_Workflow.Rmd) | Paired mouse-brain MSI/RNA, regional analysis and exported reports |

The companion-resource and extended workflows download versioned data on
first use unless a verified local or cached copy is available. Package
versions, RaMP snapshots, HGNC archives and experimental-data releases have
separate version numbers.

## Documentation and support

This repository contains the Bioconductor submission. Its vignettes document
the native workflow. The [project website](https://genomicsmachinelearning.github.io/SpaMTP/)
also includes published and developmental tutorials; historical Seurat-based
examples use a different API.

Report package problems in the [issue tracker](https://github.com/SpaMTP-project/SpaMTP/issues).
See [NEWS](NEWS.md) for changes and the
[development repository](https://github.com/GenomicsMachineLearning/SpaMTP/tree/developmental)
for website development.

## Citation

Causer, Lu, Kriel *et al.* (2026).
[SpaMTP: integrative statistical analysis and visualization of spatial
metabolomics and transcriptomics data](https://doi.org/10.1038/s41592-026-03140-8).
*Nature Methods*, 23, 1501–1506. In R, use `citation("SpaMTP")`.
