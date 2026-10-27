# Maintainer data recipes

SpaMTP loads annotation databases through SpaMTPdb. Large RaMP tables,
pathway graphs and experimental resources are not bundled in this software
package. SpaMTPdb owns annotation manifests and publication recipes;
SpaMTPData owns experimental datasets.

## Current recipes

The shipped `inst/workflows/prepare_vignette_data.R` recipe reproduces the
compact real-tissue examples in `inst/extdata/`. See `help("vignetteData",
package = "SpaMTP")` for their provenance and selection rules. The other
installed workflows prepare mouse-brain inputs and run the analysis; see
[the workflow guide](../inst/workflows/README.md).

For annotation resource updates, follow
[SpaMTPdb's maintainer guide](https://github.com/BCRL-tylu/SpaMTPdb/blob/main/CONTRIBUTING.md).
New content needs a new resource version and verified external files.

## Historical RaMP migration scripts

`update_ramp_data.R`, `update_ramp_graph_data.R`, `audit_ramp_graph_data.R`
and their helpers document the transition from older bundled RaMP objects.
They expect a populated `data/` directory containing the historical tables
and graph objects. Those inputs are absent from this submission checkout.
Run them only in an archival or disposable working copy containing the
required source objects.

The table updater reads an upstream RaMP SQLite database, retains the
previous source-identifier mapping, and uses stable identifiers to migrate
pathway graph nodes. An explicit graph-era `source_df.rda` is needed when
the graph and table snapshots differ. Graph updates support `--dry-run`;
the audit checks node membership and stable labels.

## Local output

This directory is excluded from R package builds. Its maintained scripts
remain in Git for provenance; generated `*-validation-*/` and
`native-demo-check-*/` directories are ignored. Use the ignored `.local/`
directory at the package root for other checks, rendered reports and scratch
files. Do not add large SQLite databases, downloaded resources or generated
reports to the package source.
