# SpaMTP 0.99.18

* Marker ranking is now the default in `findAllDEMs()`. Explicit technical
  pooling divides sums by the number of pixels before the log2 transform;
  unequal region sizes no longer appear as fold changes. Random pools remain
  technical subdivisions, not independent biological replicates.
* `normalizeSMData(..., "LogNormalize")` and TMM log output use log2 units,
  matching `runSpaMTPWorkflow()`. Previously saved natural-log layers require
  explicit conversion or recomputation before comparing effect sizes.
* Pre-release interface cleanup: annotation tolerances use `ppm` (default 5),
  `check_formula_mass` names the formula/ion consistency check, and
  `corr_threshold` replaces the misspelling. Functions creating an alternative
  experiment share `outputExperiment`; `outputAssay` denotes an expression
  matrix. Removed unused compatibility arguments and private helper help pages.
* `getKmeanClusters()` defaults to the `pca` reduction produced by
  `runMetabolicPCA()`. Clustering and ion plots support unnamed observations by
  position. Finite-value checks preserve sparse assays.
* Fixed unnamed numeric bin tolerances, duplicate PDF extensions, pathway
  selection without source IDs, order-sensitive named DE lists, Jaccard
  distances, and pathway thumbnails with one or two features or irregular grids.
* `workflowMarkers()`, `pathwayMapping()` and `summarizeFeatures()` expose
  results and descriptive summaries used in the guides.
* Pathway interaction preparation and historical snapshot corrections now
  belong to SpaMTPdb (>= 0.99.6), including provenance and regression tests.

# SpaMTP 0.99.17

Changes in response to the Bioconductor review of 30 October 2026.

* Annotation accuracy and usability review: CHO scoring checks formula/ion
  mass consistency instead of penalizing negative fractional mass defects.
  FMP10 product masses include the electron correction; isotope scoring
  includes multimer, solvent and derivatization atoms. Native SMILES group
  counts distinguish amines from nitro, sulfonamide and quaternary nitrogen,
  and alcohols from phosphate OH groups.
  Rebuild saved indexes to apply corrected rules. Use
  `infer_structure = "always"` to refresh older precomputed structure columns;
  `"never"` and `"always"` no longer fetch a precomputed structure table.
* Equal-mass features no longer multiply candidate rows. Index reuse honors
  adduct restrictions and rejects conflicting construction settings. Stored
  provenance comes from the reference actually used and never downloads an
  unrelated RaMP release during storage or inspection.
* Custom annotations append by default and preserve unmatched labels;
  `append = FALSE` replaces only matched values. These remain unscored display
  labels. FMP10 candidate links retain numeric reference identities, including
  masses that do not round-trip through text. Error summaries are in ppm;
  custom labels also retain `mass_diff_da`. Zero Da tolerance accepts exact
  matches and all thresholds are inclusive.
* METASPACE preserves nested isotope metadata and unavailable image URLs.
  Image conversion warns when omitting unavailable images, and rounded mass
  collisions retain unique feature names and full-precision numeric masses.
  Downloads reuse annotation image metadata, grayscale PNGs are supported,
  and hotspot clipping follows the Python client's 99th-percentile rule.
* Annotation search treats punctuation literally, refinement counts stored
  candidates when available and excludes unannotated features, and single-pixel
  QC profiles no longer fail. Numeric pathway-statistics queries require a
  match within `mz.ppm` (default 5) and reject ambiguous nearest features.
  Explicitly missing mass columns and unnamed features give actionable errors.

* Compound pathway assays follow preserved feature identities after renaming,
  subsetting or padded mass-name imports. Older stores with equal-mass rows
  and no preserved feature keys require reannotation for pathway analysis
  and candidate lookup. Annotation search uses the selected experiment's
  feature keys when display-label columns are absent.
  Normalization validates sparse
  inputs without creating a dense logical matrix.
* Removed six unused suggested packages. HTML network and density exporters
  default to the session's temporary directory; supply an output directory
  to retain the files. Analysis warnings remain visible in the gene guide,
  and bundled-resource provenance and development assistance are documented.
* Annotation parses each matched reference formula once per query and reports
  unavailable formula evidence in one warning. `FormulaStatus` identifies
  these candidates; partial element counts are excluded from formula-based
  scoring and water-loss eligibility. Scores for complete formulas are
  unchanged; affected candidates can rank higher when a penalty based on
  incomplete counts is removed. Missing evidence is recorded as `NA` and
  omitted from the score product. Target/decoy checks consolidate these warnings.
* CSV import selects dense or sparse storage from the observed density and
  preserves single-feature tables. Annotation help now states the current
  experiment defaults; installation text and vignette examples are shorter.
* The pathway viewer fits its container on resize, keeps zoom separate from
  layout, and provides a collapsible details panel. Node sizes and labels adapt
  to the view; spatial previews preserve coordinate aspect ratios. SVG export
  includes the full graph. A paired RNA/MSI map complements the live example.
* `embedPathwayNetwork()` embeds the complete offline network viewer in HTML
  documents and provides a standalone download. The pathway guides now show
  a live human RNA/MSI network, regional enrichment plots and label-free
  co-regulation; the annotation guide adds isotope and matrix-ion checks.
  Network views preserve filtered compound identities without inventing scores,
  reconcile leading-edge genes to equivalent topology records, and export
  spatial values for all coloured nodes. Missing values remain missing in
  colour scales; the viewer adapts to embedded and narrow layouts.
* Annotation candidate tables, score summaries and score plots now follow
  retained data rows after subsetting, reordering and renaming. A preserved
  feature key distinguishes equal-mass features; older ambiguous stores require
  reannotation. Custom and FMP10 annotations retain all exact nearest-mass ties
  and keep feature names, RaMP IDs, scores and stored candidates consistent.
* Vignettes have clearer introductions, corrected cross-references, formatted
  result tables and a shared responsive layout. Memory guidance distinguishes
  delayed scaling from dense allocations during graph PCA and modality merging.
* Spatial feature plots accept unnamed pixels and preserve column alignment.
* Enrichment, SMILES parsing and METASPACE image downloads use `BPPARAM`.
  Serial execution is the default; progress follows `verbose` unless an
  explicit parameter object is supplied. This replaces `nproc`, `workers`
  and `structure_workers` and never changes the session's future plan.
* The variance plot requested by `runMetabolicPCA(show_variance_plot = TRUE)`
  is stored in `metadata()$reduction_plots` instead of printed.
* `densityMap()` exports an offline viewer using the plotly R package and
  bundled d3. Copy its HTML and `density_map_files` directory together.
  d3's copyright and ISC license are recorded in `inst/COPYRIGHTS` and
  `inst/vendor/d3/LICENSE`; only the minified distribution is shipped.
* Formula masses and element counts use MetaboCoreUtils; adduct groups use
  igraph connected components. Removed the unused legacy `adduct_file`.
  Annotation edges are assembled in batches, SMILES parser storage is
  preallocated, and workflow callbacks use private record collectors.
* Examples now execute small offline analyses using `spamtp_example`, with
  explicit interactive guards for the remote METASPACE service. The shipped
  GitHub installer and its installation test have been removed.
* Removed unused or replaceable dependencies (`purrr`, `mclust`, `future`,
  `future.apply`, `parallel`, `magick`, `tibble`, `naturalsort`). DropletUtils,
  scater, pheatmap and ape are optional and checked where they are needed.
* Pathway coverage matches memberships in batches; m/z queries locate all
  mass windows in one pass. The annotation vignette queries a deterministic
  feature subset against the complete reference, keeping full-panel MS1
  context and normalization totals.

* Breaking: arguments that name an expression assay are now called
  `assayName` (previously `slot`), with `smAssayName`, `stAssayName`,
  `mzAssayName` and `pathwayAssayName` replacing `SM_slot`/`SM.slot`,
  `ST_slot`/`ST.slot`, `mz.slot` and `pathway.slot`, and `assayNames`
  replacing `slots` in `plot3DFeature()`. `convertSeuratToCardinal()` takes
  `layer`, matching the Seurat 5 term.
* New plotting and inspection functions, so analyses no longer need
  hand-written plotting or helper code:
  - `plotSpatialLabels()` draws categorical pixel annotations at their spatial
    positions.
  - `plotSpaMTPWorkflow()` returns the figures of a `runSpaMTPWorkflow()`
    result as ggplot objects.
  - `annotationCandidates()`, `summarizeAnnotationScores()` and
    `plotAnnotationScores()` rank, summarize and plot stored annotation
    candidates.
  - `markerRanks()` arranges `findAllDEMs()` statistics as a feature-by-region
    rank matrix.
* `plotSpatialFeature()` has a readable default theme and a `legendTitle`
  argument, and finds features in paired alternative experiments, so ion and
  transcript maps can be drawn together.
* `plotPathwaysSpatially(combine = TRUE)` arranges pathway maps with a shared
  colour bar.
* `findNearestMZ()` accepts several target masses (`mz`, previously
  `target_mz`).
* `findRegionalPathways()` gains `seed`, applied with `withr::local_seed()` so
  the caller's random number stream is unchanged, and `regions`. Its `ranks`
  argument accepts `findAllDEMs()` results, and the result reports measured
  and leading-edge genes and metabolites per pathway.
* `pathwayNetworkPlots()` attaches a `network_summary` table to its result.
* PCA uses `scrapper::runPca()` in `runMetabolicPCA()` and
  `multiOmicIntegration()`, and UMAP displays use `scrapper::runUmap()`.
  BiocSingular is no longer imported. The 500 most variable features are
  still used by default (`ntop`).
* Matrix storage:
  - `normalizeSMData()` stores mostly nonzero results as an ordinary
    `matrix` rather than a dense-valued `dgCMatrix`, keeps sparse results
    sparse, and keeps `DelayedArray` input delayed.
  - Its new `libSizes` argument normalizes by recorded totals, such as the
    full-spectrum TIC of a feature subset.
  - `scaleSMData()` returns a lazily evaluated `DelayedMatrix` by default,
    using block-processed row statistics, and documents the memory needed for
    a materialised result.
* Adduct ion masses come from `MetaboCoreUtils::adducts()`. SpaMTP adds its
  rule priors and structural constraints, and records the source definition in
  `adductRules()$metabocoreutils_name`.
* Removed `utils::globalVariables()`: data-frame columns in ggplot2 and dplyr
  code use the `.data` pronoun or quoted names, so no workspace variable can
  be picked up by mistake. NAMESPACE no longer imports whole packages; every
  import is selective and no symbol is imported twice.
* Reworked result-building loops in `findRegionalPathways()`, pathway networks,
  Fisher input mapping and the workflow report to combine intermediate records
  once. Two unused internal functions were removed, so MSnbase is no longer
  imported.
* `getMetaspace()` and `refineLipids()` gain `verbose`; enrichment progress
  follows the default `BPPARAM` settings described above.
* `bladder_demo.rds` includes the archived `compartment` labels.
* Vignettes:
  - They call package functions instead of building plots or helper functions.
  - They attach packages instead of using `pkg::` calls.
  - They contain no unit checks, loops or global seeds.
  - All figures are generated at run time.
  - Installation instructions use Bioconductor only.

# SpaMTP 0.99.16

* Shortened the README and clarified the introduction, staging installation,
  companion package roles and reading order of the seven vignettes.
* Moved detailed installation guidance into `inst/workflows/README.md` and
  excluded local development artifacts from package builds.
* Added opt-in annotation quality control. Defaults of `annotateSM()` and the
  index functions are unchanged.
  - `estimateAnnotationFDR()` scores decoy copies of a peak list moved by a
    constant mass, which keeps isotope and adduct spacings exact, and reports
    a decoy-to-target match ratio at each score cutoff;
    `annotationScoreCutoff()` picks a cutoff for that diagnostic ratio.
    These diagnostics do not establish calibrated false-discovery control.
  - `flagIsotopePeaks()` marks 13C isotopologues by spacing, intensity ratio
    and co-localisation, and `flagMatrixPeaks()` / `matrixIonTable()` mark
    DHB, CHCA, 9-AA, DAN and norharmane monomer and cluster ions.
  - `maldiAdductPreset()` returns a matrix rule set without multiply charged
    and electrospray-solvent adducts.
  - `groupAdductPeaks()` links co-localised peaks at exact adduct spacings
    (and split peaks), optionally excluding flagged isotope and matrix peaks,
    and `annotateAdductGroups()` gives each grouped peak the compounds that
    explain the most of its group's ions through at least two adducts.
  - `annotationDatabaseScope()` restricts a RaMP chemical table to compounds
    in selected pathway sources.

# SpaMTP 0.99.15

* Rendered vignette figures at 96 dpi (192 dpi on high-density displays) and
  excluded the README overview figures from the source tarball, bringing the
  built package under the Bioconductor 10 MB limit.
* Installer tests no longer depend on the R version running them; a separate
  test covers the R >= 4.6.0 guard.
* Regional marker groups, workflow marker tables and association ties are now
  ordered the same way in every locale. `R CMD build` sorts with
  `LC_COLLATE=C`, which previously changed the workflow's association anchors
  and failed the end-to-end vignette.

# SpaMTP 0.99.14

* Streamlined vignette interpretation while retaining observed metrics,
  pathway significance and identity ambiguity. Removed an artificial
  inconsistent-mass fixture from the user guide; regression tests retain
  coverage of the mass-defect penalty.
* Reconciled vignette descriptions with executed code and displayed output:
  equal-weight integration, visible QC summaries, pinned native-data retrieval,
  Fisher table columns, pathway names and topology repair provenance.
* Added positive single- and double-water-loss rules gated by alcohol groups,
  neighbouring hydrogen availability and valence-compatible losses. Older
  structural tables are checked only for mass-matched candidates; each isomer
  retains its own eligibility. Candidate tables expose the capacity and
  eligibility status.
* Demonstrated engine-generated sphingosine dehydration candidates within the
  complete RaMP search and a formula-specific structural comparison. Clarified
  annotation ambiguity across score cutoffs and the difloxacin mass candidate.
* Explained the human example's Xenium observation unit, log transformation
  bases and joint-pathway coverage. Standardized terminology, abbreviations,
  section links and numeric formatting across the seven guides.

# SpaMTP 0.99.13

* Numbered all seven vignette titles and index entries to guide readers from
  the introductory workflow through data access, companion resources,
  metabolite annotation, gene mapping and pathway databases to the complete
  spatial multi-omics workflow.
* Added internal subsection navigation, declarative explanations, contextual
  parameter definitions and direct prose. Refined table spacing and figure
  layout while retaining BiocStyle typography and colours. Credited Andrew Causer and
  Tianyao Lu in all vignette headers.
* Numbered figures by vignette (for example, Figure 1.1 and Figure 7.1) with
  linked in-text references. Added concise guidance on input alignment, stored
  provenance and interpretation of annotation and pathway results.
* Added figure captions and consistent display sizing. The complete workflow
  now presents a compact regional pathway plot with wrapped labels and separate
  legends, while retaining the complete pathway results and adjusted P values.
  These layout defaults are shared by `plotRegionalPathways()`, with an optional
  `label_width` argument for controlling name wrapping.
* Replaced tiny spatial demonstrations with tissue-scale examples. The bundled
  bladder panel retains 2,210 positions on a four-pixel grid, archived SSC labels
  and explicitly historical lipid candidates. A new pig fetal DESI panel keeps
  all 4,959 locations and three prespecified ion windows, normalized with the
  full 10,200-bin total ion current. Source checksums, licenses and preparation
  recipes are recorded; archived labels and candidates are not treated as
  independent histological truth or confirmed chemical identifications.
* Mouse-brain vignettes use all 876 paired DHB/Visium locations from pinned
  SpaMTPData 1.1.0 to preserve lesion and anatomical context. Removed the unused
  160-location bundled DHB subset (1.45 MB).
  Report headings and figure titles use declarative wording.

* Corrected annotation headings, chemical structure rendering, pathway links and
  mass-defect result reporting. Added score distributions, candidate identities,
  ambiguity counts at the downstream threshold, and network size summaries.
* Normalized RNA before pathway maps and used regional markers for Fisher
  enrichment. Added a human MS+RNA workflow with modality coverage and explicit
  interpretation of its limited compound support and nonsignificant results.
* Documented the human glioma study, acquisition platforms and source archives.
  Distinguished archived processed MS values from RNA counts; bladder plots
  retain upstream TIC scaling and include a panel-normalization sensitivity check.
* Interpreted clustering agreement, modality silhouettes and covariate-adjusted
  associations. Ranked displayed markers by AUC and qualified the sphingosine
  mass match as a hypothesis for chemical follow-up.
* Executed a temporary Cardinal imzML round trip, added installation and citation
  guidance, documented function-name migration, and made pathway topology
  provenance available through installed resources and a compact figure.

# SpaMTP 0.99.12

* Added a standalone installation workflow with a longer download timeout,
  explicit remotes setup, official SpaMTPData installation and a staging
  fallback for SpaMTPdb while it awaits repository publication.
* Vignette installation checks Pandoc before installing the full dependency
  set and discovers common Windows RStudio/Pandoc locations. Installation
  settings propagate to child R processes and are restored on completion or
  failure. Windows troubleshooting distinguishes incomplete ZIP downloads,
  unavailable repositories and missing document tools.
* Report export removes a verified empty destination before publishing the
  completed report, accommodating Windows directory-rename restrictions.
  Files added to the destination during rendering are preserved.

# SpaMTP 0.99.11

* Vignettes use bundled, reproducibly selected real spatial subsets: 160
  DHB/Visium locations and 240 human-brain locations. All measured features,
  raw counts, paired modalities, coordinates and source labels are retained.
  Full-section mouse-brain analyses remain available as installed recipes.
  All vignette analysis chunks execute during builds.
* Feature correlations read sparse assays in bounded blocks and residualize
  multiple features together, retaining raw and adjusted Pearson results and
  undefined-feature handling.
* Graph PCA uses matrix-vector products for wide inputs requesting few
  components, avoiding a dense feature-by-feature covariance matrix while
  retaining the same graph-regularized objective.

# SpaMTP 0.99.10

* Added return-value documentation to the experimentAccess guide, resolving
  the missing-value-section warning from BiocCheck.

# SpaMTP 0.99.9

* Clarified native data access across function help pages: modalities use
  altExp(), expression matrices use assay(), and retained slot arguments
  select matrix names rather than S4 storage slots. Corrected return types,
  reduction locations and normalization descriptions. Added a runnable
  experimentAccess guide. Assay selectors now require exact character names;
  modality merging rejects duplicate primary aliases and prefers native
  scaled values over the legacy scale.data matrix.

* All vignette R chunks now execute by default, including full published
  human-brain, FMP10/Visium and DHB/Visium analyses. Resource versions are pinned;
  downloads and long builds are expected when verified caches are absent.
  Original Space Ranger scale factors supply the registered FMP10 pair's spot
  geometry. Synthetic file-I/O and replicate-design examples remain labelled.
  Installed mouse-brain recipes can resolve public resources when no local
  directory is supplied, while explicit directories retain offline defaults.

* Annotation statistics now select the requested MSI modality's candidate store
  for both single-feature and batch ranking, without borrowing another modality's
  latest root-level annotations. Historical and assay-local stores remain usable.
* Regional marker scores, spatial-block sensitivity and biological-replicate
  contrasts preserve matrix dimensions for single-feature inputs.

* Replaced the installed mouse-brain demo and preparation recipes with native
  SpatialExperiment workflows. Installed examples no longer load a legacy
  checkout, read S4 slots directly or require Seurat. Compound networks reuse
  scored annotations, scran-based regional effects and a shared pathway index.
  Optional optical images and spot geometry are supplied explicitly; FMP10
  mapping refuses a missing radius rather than guessing it from pixel spacing.
  Existing output paths are preserved and script sourcing has no side effects.

* Fixed pathway interaction corruption in the graphite-derived RaMP 3.0.7
  graphs. A pinned, checksum-guarded correction restores per-edge source labels,
  directions, and parallel interactions for KEGG, Reactome, WikiPathways and
  SMPDB. Interaction styles now use label semantics across databases, and
  undirected edges have no arrowhead. The original published resource bytes
  retain their checksums; SpaMTP applies corrections when loading them.
* Added the Pathway Database Integration vignette, a reproducible correction
  builder, source checksums, and regression tests for the Cell Cycle report.


* Reorganized the workflow around measurement quality, correspondence,
  representation comparisons, regional markers, replicate inference,
  conditional cross-omic association and shared-identity pathway analysis.
  Reports connect native heatmaps, H&E feature maps, effect plots, 3D views
  and pathway networks with their actual analytical inputs.
* Native PCA accepts explicit features and records the fitted feature set.
  Graph PCA uses a sparse Laplacian and all retained observations; its
  embeddings now feed spatial joint integration and optional clustering.
  Comparisons use matched feature selection/scaling, multiple requested K,
  external-reference ARI, approximate silhouettes and block-omission
  clustering sensitivity conditional on the fitted embeddings.
* findAllDEMs() supports descriptive scran marker effects/AUC and explicit
  biological-replicate limma contrasts. Spatial-block marker sensitivity
  carries no invented P values. The historical technical-pool mode warns
  that its pools are not biological replicates. Native heatmaps preserve
  the actual analysis matrix and rank numeric FDR correctly.
* findCorrelatedFeatures() supports generic feature IDs, explicit target
  sets and raw versus covariate-residual correlations. Exact RNA/protein
  features no longer require m/z metadata for native plots.
* findRegionalPathways() accepts complete unfiltered identity-level ranks.
  Competitive enrichment, scores, GESECA and network views reuse the same
  membership index. Compound identity assays support sparse inputs and
  audit ambiguity exclusions and exact aggregation weights. Annotation
  stores are preserved separately for multiple MS modalities.
* plotRegionalPathways() displays BH-adjusted significance and orders
  pathways using 1 minus Jaccard similarity, including single-pathway input.
  The region explorer adds directional AUC and marker-stability filters.
* Reports audit pathways with identical measured members across regions without
  shrinking the BH family. Native networks embed a pinned, licensed D3
  distribution for offline use, retain missing effects as unavailable, and
  label descriptive effects in the declared workflow units.
* Added portable FMP10/Visium and DHB/Visium recipes using pinned
  SpaMTPData/SpaMTPdb resources, plus a rewritten analytical tutorial with
  three synthetic modalities. Mouse genes are not mapped to human HGNC
  by changing case. Equal-weight joint embeddings remain distinct from
  the historical Seurat WNN method.

# SpaMTP 0.99.8

* Added analyzeSpaMTPRegions() and a portable offline region/DE explorer.
  Modality, region and replicate contrast selectors link feature tables,
  effect/volcano plots, spatial maps, regional means and distributions.
  Search, effect/direction/FDR filters, point inspection, zoom/pan and filtered
  CSV export operate on an explicitly labelled bounded preview; full results
  are exported separately. Blank region labels are excluded and audited.
* Reports now show native plotSpatialFeature(), mzViolinPlot(),
  findSpatiallyVariableMetabolites() and runSpatialGraphPCA() outputs, plus
  regional pathway score heatmaps when a compatible index is configured.
  Native spatial analyses record selected features and positions per sample.
  Descriptive region effects never receive pixel-based DE P values.
* Method descriptions use smaller serif text and collapsible detail panels.
  Capability cards identify computed and unavailable native modules.
  Existing workflow results can be extended without changing their original
  contrasts; core, extension and renderer versions remain distinguishable.
* Expanded the generic workflow tutorial with paired mouse brain, independent
  single-omic inputs, arbitrary region metadata, preview controls and native
  analysis settings. No study-specific feature names or species conversion
  are embedded in the workflow or report renderer.
* Single-modality workflows now apply requested clustering to their PCA and
  expose those exploratory labels to the region browser. Previously, clusters
  were only computed after joint integration and a single-omic request was
  silently ignored. Cluster input and settings are recorded.

# SpaMTP 0.99.7

* Added runSpaMTPWorkflow() for explicit paired or independent spatial
  modalities: pinned companion data acquisition, QC, registration/mapping,
  normalization, per-modality PCA, joint representation, exploratory
  associations, biological-replicate contrasts and optional annotation/pathways.
  Unmatched positions are audited and excluded from joint analysis.
* Added renderSpaMTPReport() with embedded figures, full CSV tables, native
  result RDS, registration weights, resource provenance and file checksums.
  Rendering needs neither Pandoc nor network access and preserves existing files.
* Added an offline end-to-end vignette demonstrating two/three modalities,
  known-transform registration and portable reporting, with explicit boundaries
  for raw preprocessing, biological inference and platform-specific settings.

# SpaMTP 0.99.6

* Added buildPathwayIndex() to share reconciled gene identities and pathway
  memberships across enrichment, RaMP expression assays, pathway scores,
  GESECA, named plots and network coverage. Pathway IDs distinguish different
  sets with the same display name. Reuse validates membership fingerprints.
* createPathwayAssay() now maps human genes with the audited HGNC index and
  preserves the selected expression layer. Identical duplicate gene rows are
  counted once; different rows require explicit mean or sum aggregation.
  createPathwayObject() accepts original gene identifiers or mapped assays.
* Pathway outputs retain database sizes before/after reconciliation, measured
  coverage, actual used members and excluded conflicting RaMP records.
  Zero-coverage and size-excluded pathways remain visible in coverage audits.
  Scoring and enrichment use the same measured membership sets. Regional and
  GESECA size filters use members actually supplied to their test engines.

# SpaMTP 0.99.5

* Added `buildGeneMappingIndex()`, `mapGeneIdentifiers()` and native
  `annotateGeneIdentifiers()` using SpaMTPdb's versioned HGNC reference.
  Approved/previous/alias symbols and stable IDs have explicit resolution
  status; ambiguous aliases and conflicting RaMP records are never expanded
  into multiple genes. Gene-specific IDs can resolve shared protein IDs.
* Fisher and regional pathway analyses merge RaMP records belonging to the
  same HGNC gene and unite their pathway memberships before counting. Network
  displays reuse the same identity mapping. Conflicting differential values
  from multiple features of one gene require resolution before analysis.
  Gene mapping and reference provenance are retained in outputs and rowData.
* Gene indices now record a normalized source-table fingerprint and reject
  incompatible reuse even when RaMP record IDs are unchanged. Older index
  schemas must be rebuilt. Shared crossreferences retain unresolved and
  conflicting raw candidates instead of implying a unique match. Conflict
  membership is indexed once and repeated queries are resolved once while
  preserving the original input order and per-row audit.
* Official database workflows use HGNC mapping automatically for human genes;
  custom fixtures remain offline, and `gene_mapping = "ramp"` supports explicit
  historical reproduction or curated non-human resources. HGNC mapping does
  not perform cross-species orthology conversion.

* Fixed the Fisher pathway contingency table: foreground non-members are
  K - overlap, not max(0, K - pathway_size). Added an explicit measured
  `universe`, deduplicated RaMP memberships and validated foreground inclusion.
  Background size is fixed before pathway-size filtering. Every eligible
  pathway, including zero-overlap sets, now belongs to the BH testing family.
  `pval_cutoff` now filters FDR as documented (previously it filtered raw p).
  Results include foreground/background counts, internal pathway IDs and an
  `enrichment` audit attribute; missing display metadata no longer drops tests.
  These corrections intentionally change previous p-values and FDRs.

# SpaMTP 0.99.4

* Updated the companion workflow for SpaMTPData >= 0.99.4 and its published
  native resource release 1.1.0. Historical archive examples explicitly select
  1.0.0; native downloads require no Seurat conversion.
* Preserve non-syntactic feature and pixel metadata names during Seurat
  conversion, including labelled metabolite intensities in the human brain
  archive. Added forward/reverse conversion regression tests.
* Support the companion data package's one-time native resource preparation:
  resulting SpatialExperiment RDS files are used directly without Seurat.

# SpaMTP 0.99.3

* Synchronized companion-package calls with the camelCase SpaMTPdb >= 0.99.2
  and SpaMTPData >= 0.99.1 APIs. Added an evaluated three-package workflow
  using a native SpatialExperiment example without downloads or Seurat.
* Passed local-file checksum verification through loadSpaMTPDatabase(), and
  isolated in-session cache entries by configured resource directory and
  verification setting. Explicit custom database bundles remain supported.
* Updated the full annotation vignette to use the current registry's bytes
  field rather than the obsolete serialized_bytes field.

# SpaMTP 0.99.2

* Migrated the complete indexed annotation tutorial to native containers and
  added a runnable public-accessor guide for expression, metadata, coordinates,
  images and paired modalities. Static regression checks cover actual R code
  in both functions and vignettes, including unevaluated examples, without
  mistaking SMILES chirality or `slot=` argument names for slot operations.
* Seurat spatial conversion now prioritizes a single image/FOV accessor over
  historical metadata coordinates and records the selected coordinate source.
  Multiple images require explicit selection; finite-radius centroids are read
  as pixel centres rather than expanded polygon vertices. Explicitly supplied
  coordinates remain supported, and metadata is the fallback for image-free
  objects only. Input objects and their archival metadata are preserved.
* Native CSV and METASPACE imports store coordinates only in `spatialCoords()`,
  removing redundant metadata copies that could become stale after alignment.
* Made `SpatialExperiment` the default container for aligned/binned MSI and
  retained Cardinal `MSImagingArrays` / `MSImagingExperiment` for raw, continuous, or file-backed
  spectra. Added S4 conversion generics, registered Cardinal/SPE coercions,
  and explicit optional Seurat conversion functions.
* Added raw-array S4 preprocessing methods, execute deferred Cardinal
  processing before binning/conversion, and cover the mass range of all raw
  spectra rather than inferring it only from the first pixel.
* Mapped MSI feature annotations, pixel metadata, coordinates, optical images,
  and paired transcriptomes to `rowData()`, `colData()`, `spatialCoords()`,
  `imgData()`, and `altExp()`, respectively. Seurat was removed from dependency
  declarations; SeuratObject is optional and used only by explicit converters.
* Added S4 methods for binning, normalization, multi-omic integration, spatial
  plotting, images, and paired transcriptomes. Native workflows reuse
  Cardinal, SpatialExperiment, SingleCellExperiment, scater, edgeR, limma,
  DropletUtils and fgsea functionality.
* Removed unused direct Imports of EBImage, sp, shinyjs, RColorBrewer, zeallot
  and matter after migrating the image/alignment workflows. Raw MSI continues
  to use matter through Cardinal; optical images use SpatialExperiment.
* Pathway and merged-modality matrices are stored as alternative experiments,
  allowing their feature dimensions to differ from the primary MSI assay.
* Migrated affine alignment, image attachment, pixel/spot mapping, ROI selection,
  spatial correlations, Moran's I and pathway-score plots to native containers.
  Mapping and spatial graphs separate sample identities; graphs also use
  `colPairs()` so standard subsetting updates their indices.
* Added regression coverage for reordered pixels, multi-sample coordinates,
  altExp normalization, intensity-conserving matrix bins, image scale factors,
  constant features and RaMP feature IDs. Fixed limma treatment-table extraction
  and preserved explicit sample labels during technical pooling with edgeR.
* Native integration defaults to main MSI plus the paired transcriptome; derived
  pathway and merged assays are not integrated automatically. Equal-weight
  concatenated PCA is explicitly distinguished from historical Seurat WNN.
* Migrated MSI feature/annotation plots, optical overlays, 3D views, mass
  spectra, interactive mass windows, density export and pathway-network
  container extraction to SpatialExperiment/altExp and imgData.
  Plotting rejects unsupported Seurat-only arguments instead of ignoring them.
  Mass windows now sum all in-range features and overlapping windows count
  each peak once; m/z axes no longer require specially formatted feature IDs.
  Fixed row-major raster colour ordering for transparent optical overlays and
  3D image points; S4 container classes are now imported explicitly for coercion.
* Replaced pathway-based annotation ranking's temporary Seurat/Cardinal
  conversion with direct paired Pearson correlations. Reused existing pathway
  scoring and current RaMP candidate resolution; respected requested assays,
  removed duplicate candidate IDs, and defined constant/missing-data handling.
  The retained weighted z-score tail probabilities are explicitly documented
  as heuristic ranks, not calibrated identification p-values.
* Fixed annotation candidate sorting when optional Score/Error columns are
  absent. Retired four unused internal plotting helpers and their help pages.
* Removed non-conversion Seurat branches and implicit Seurat importer output.
  Added `seuratToSingleCellExperiment()` for non-spatial input, explicit named
  coordinate conversion, strict Seurat layer matching and identity transfer.
* Added `scaleSMData()` for feature-wise centring/scaling, and corrected PCA
  and expression plotting to use the requested alternative experiment.
  QC plot inputs and `runDE()` now use `data`, not Seurat-named parameters;
  `subsetSPM()` no longer accepts retired object/slot-upgrade arguments.
* Native annotation, curated FMP10 annotation and feature subsetting preserve
  arbitrary feature IDs and paired modalities. Pathway entry points use
  native assays instead of Seurat-style double-bracket indexing.
* Added a fresh-process native workflow test that requires neither Seurat
  namespace, plus adapter-boundary and conversion regression tests.

# SpaMTP 0.99.1

* Standardised all exported function names to lower camel case for the
  Bioconductor API. This is an intentional breaking change on the submission
  branch; the published-work development branch retains its historical API.
* Removed direct S4 slot access from package and test code. Internal container
  access now uses public SeuratObject, SummarizedExperiment, S4Vectors,
  Cardinal, and MSnbase accessors. Annotation and alignment provenance is
  stored with `SeuratObject::Misc()` while legacy `Tool()` entries remain
  readable.
* Added an internal container-access layer for Seurat and
  SummarizedExperiment-derived objects, including SingleCellExperiment and
  SpatialExperiment, as the foundation for broader Bioconductor container
  interoperability.
* Removed the unsupported manual mutation of Seurat modality-weight internals
  from `multiOmicIntegration()`; weighted-nearest-neighbour construction now
  delegates entirely to `Seurat::FindMultiModalNeighbors()`.

# SpaMTP 0.99.0

* Split versioned resources from the analysis code: `SpaMTPdb` now supplies
  the pruned RaMP annotation and pathway snapshot, and `SpaMTPData` provides
  named access to large experiment/vignette objects. `loadSpaMTPDatabase()`
  and `spaMTPDatabaseInfo()` form the public database interface. The software
  package retains only the small adduct and reaction-style constant tables.
* Updated the package citation to the peer-reviewed *Nature Methods* paper,
  designated Tianyao Lu as the current package maintainer, and added direct
  links to the developmental documentation and source branch. Andrew Causer
  remains credited as an original author and former maintainer.
* Added `applySpatialAlignment()` to apply existing SMINT coordinate outputs or
  homogeneous transforms, estimate landmark-based affine alignment in R, and
  run the SMINT-compatible STalign LDDMM workflow through an optional Python
  backend. Alignment provenance and nearest-target diagnostics are retained in
  the returned SpaMTP object without storing large Python velocity tensors.
* Pathway enrichment, pathway-assay construction, Fisher analysis of m/z
  inputs, and interactive pathway networks now consume the scored `Ramp_IDs`
  produced by the indexed annotation engine. `annotateSM()` resolves the
  versioned RaMP 3.0.7 chemical-property table through `SpaMTPdb`, records
  annotation provenance through `SeuratObject::Misc()`, with a compatibility
  fallback for legacy serialized objects.
  Legacy annotations require an explicit `annotation_source` fallback.
* Annotation candidates can now be stored once with `annotateSM(min_score =
  0)` and filtered later with `annotation_score_threshold`. Interactive
  pathway networks add `metabolite_detection = "annotated"` to display
  score-filtered pathway metabolites without requiring DE significance, while
  retaining separate leading-edge evidence in node styling and tooltips.
* Moved the pruned RaMP snapshot from 2.5.4 to versioned SpaMTPdb 3.0.7
  resources, with explicit upstream/source metadata and reproducible resource
  staging scripts.
* Added a chemically validated adduct rule table and reusable sorted m/z index.
* Added ppm pruning, proton/charge bounds, mass-defect scoring, isotope checks,
  and contextual adduct-family scoring while retaining the legacy annotation
  result columns.
* Added dependency-free SMILES graph decomposition, functional-group and atom
  site reporting, mode-specific protonation/deprotonation and alkali-binding
  priors, and `predictAdductsFromSMILES()`. Aromatic functional groups are
  perceived consistently from both lowercase aromatic and alternating-bond
  six-membered Kekule SMILES. Full RaMP runs join an independent precomputed
  `SpaMTPdb::smiles_features` resource; small custom databases are inferred
  automatically at runtime.
* `chem_props` can now be used directly as a RaMP-backed annotation database.
* Rebuilt `pathwayNetworkPlots()` around cached topology lookup, precomputed
  edges, selective sparse-matrix extraction, and JSON serialization. The new
  responsive viewer adds zoom/pan, focused/full networks, label controls,
  shared reaction markers, spatial inspection, and SVG export.
* Made legacy pathway annotation robust to unequal `Isomers_IDs` and
  `IsomerNames` counts, and added a cached Mouse Brain/DHB network demo.
* Pathway network metabolite labels now remain visible in the default label
  mode; dragged nodes stay pinned until released, and edge legends persist
  when the focused view changes.
* Added interactive repulsion, balanced-force, radial, and gene/metabolite
  pathway layouts, with a more widely spaced repulsion layout as the default.
* Harmonised all versioned pathway-graph node IDs with RaMP-DB 3.0.7, added a
  reproducible cross-version graph updater/auditor, and added stable source-ID
  labels when RaMP itself has no common name.

# SpaMTP 1.1.0

***Updated SpaMTP Release (Oct 2025)***

Additional Features:

* Updated `loadSM()` function for compatability with **Cardinal V3.8**
* Implementation of GraphPCA in R.
* Additional functions for handling large datasets - `annotateBigData` and `selectROIs`.
* Function for refining m/z annotations based on correlated pathway activity or MS/MS profiles.
* Package wide update for compatability with **Cardinal V3.8** and **Seurat V5.3**.


# SpaMTP 1.0.0

* Initial SpaMTP Release (March 2025)
