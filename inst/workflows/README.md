# SpaMTP workflow recipes

Install SpaMTP and its example data from Bioconductor:

```r
BiocManager::install(c("SpaMTP", "SpaMTPData"))
```

## Installed workflow recipes

| File | Purpose |
|---|---|
| `prepare_mouse_brain.R` | Prepare native mouse-brain inputs with explicit spatial configuration. |
| `run_mouse_brain.R` | Run the analysis and produce an HTML report with exported tables. |
| `prepare_vignette_data.R` | Reproduce the compact measured datasets used by the introduction. |

Inspect each script's usage before running it. The numbered SpaMTP vignettes
explain the analysis choices and resource versions. Find these recipes after
installation with `system.file("workflows", package = "SpaMTP")`.
