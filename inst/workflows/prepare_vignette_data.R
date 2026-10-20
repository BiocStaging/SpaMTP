# Recreate the bundled real-data demonstrators from verified SpaMTPData 1.1.0.
# Sourcing only defines functions; call prepareVignetteData(outputDir, nativeDir).
# Selection uses coordinates and archived labels, never expression or results.

prepareVignetteData <- function(outputDir, nativeDir = NULL,
                                offline = !is.null(nativeDir)) {
    if (length(outputDir) != 1L || is.na(outputDir) || !nzchar(outputDir) ||
        file.exists(outputDir)) stop("Supply a new output directory.")
    dir.create(outputDir, recursive = TRUE)
    read <- function(resource) SpaMTPData::spaMTPData(resource,
        version = "1.1.0", local_dir = nativeDir, offline = offline)
    subset <- function(object, keep, fields, rule, parameters) {
        rna <- SingleCellExperiment::altExp(object, "SPT")
        result <- SpatialExperiment::SpatialExperiment(
            assays = list(counts = SummarizedExperiment::assay(
                object, "counts")[, keep, drop = FALSE]),
            rowData = S4Vectors::DataFrame(
                mz = SummarizedExperiment::rowData(object)$mz,
                row.names = rownames(object)),
            colData = SummarizedExperiment::colData(object)[keep, fields,
                drop = FALSE],
            spatialCoords = SpatialExperiment::spatialCoords(object)[keep, ,
                drop = FALSE])
        paired <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = SummarizedExperiment::assay(
                rna, "counts")[, keep, drop = FALSE]),
            rowData = S4Vectors::DataFrame(row.names = rownames(rna)))
        provenance <- S4Vectors::metadata(object)$SpaMTPData
        S4Vectors::metadata(result)$SpaMTPData <- provenance
        S4Vectors::metadata(paired)$SpaMTPData <-
            S4Vectors::metadata(rna)$SpaMTPData
        SingleCellExperiment::altExp(result, "SPT") <- paired
        S4Vectors::metadata(result)$vignette_subset <- list(
            source = provenance, original_observations = ncol(object),
            rule = rule, parameters = parameters, barcodes = colnames(result),
            features = c(metabolomics = nrow(object), transcriptomics = nrow(rna)),
            retained = "All measured features, archived assay values, pairing and original coordinates",
            omitted = "Archived annotations, normalized layers, reductions and images",
            interpretation = "A spatial subset of one specimen for descriptive analysis")
        result
    }
    human <- read("human_brain_integrated")
    xy <- SpatialExperiment::spatialCoords(human)
    center <- apply(xy, 2, stats::median)
    distance <- rowSums(sweep(xy, 2, center)^2)
    keep <- sort(order(distance, colnames(human))[seq_len(240L)])
    result <- subset(human, keep, c("sample_id", "cout_labels"),
        paste("Take the 240 locations nearest the coordinate-wise median;",
              "break distance ties by barcode and retain original column order."),
        list(observations = 240L, center = center))
    S4Vectors::metadata(result)$study_provenance <- list(
        tissue = "Low-grade glioma; one spatial subset of an aligned section",
        study = "https://doi.org/10.1038/s41591-025-03884-4",
        paired_data = "https://doi.org/10.6019/S-BIAD1426",
        processed_archive = "https://doi.org/10.5281/zenodo.18282678",
        metabolomics = "Positive-mode MALDI, iMScope QT, CHCA matrix",
        transcriptomics = "Xenium targeted 339-gene panel; aligned serial section",
        assay_values = paste("MS preserves processed, spatially aligned abundances;",
            "RNA preserves integer counts from the source resource"),
        cout_labels = paste("Archived cell-state groups; the exact assignment",
            "procedure is unavailable in the resource"))
    saveRDS(result, file.path(outputDir, "human_demo.rds"), compress = "xz")
    prepareBladderVignetteData(file.path(outputDir, "bladder_demo.rds"),
        nativeDir = nativeDir, offline = offline)
    invisible(outputDir)
}

# A regular coordinate grid retains the bladder section's spatial extent.
# All 79 features in the archived annotated panel are retained; this is not
# the full raw imzML feature space. No annotation or expression guides sampling.
prepareBladderVignetteData <- function(outputFile, nativeDir = NULL,
                                       offline = !is.null(nativeDir)) {
    if (length(outputFile) != 1L || is.na(outputFile) || !nzchar(outputFile) ||
        file.exists(outputFile)) stop("Supply a new output file.")
    object <- SpaMTPData::spaMTPData("mouse_bladder_annotated", version = "1.1.0",
        local_dir = nativeDir, offline = offline)
    xy <- SpatialExperiment::spatialCoords(object)
    stopifnot(all(c("x", "y") %in% colnames(xy)), all(is.finite(xy)),
              all(xy == round(xy)))
    origin <- apply(xy[, c("x", "y"), drop = FALSE], 2, min)
    keep <- which((xy[, "x"] - origin["x"]) %% 4 == 0 &
                  (xy[, "y"] - origin["y"]) %% 4 == 0)
    stopifnot(length(keep) > 0L)
    result <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = SummarizedExperiment::assay(
            object, "counts")[, keep, drop = FALSE]),
        rowData = S4Vectors::DataFrame(
            mz = SummarizedExperiment::rowData(object)$mz,
            archived_lipid_candidates =
                SummarizedExperiment::rowData(object)$Species.Name.Simple,
            row.names = rownames(object)),
        colData = SummarizedExperiment::colData(object)[keep, c("sample_id", "ssc"),
            drop = FALSE],
        spatialCoords = xy[keep, , drop = FALSE])
    provenance <- S4Vectors::metadata(object)$SpaMTPData
    S4Vectors::metadata(result)$SpaMTPData <- provenance
    S4Vectors::metadata(result)$vignette_subset <- list(
        source = provenance, original_observations = ncol(object),
        rule = paste("Keep pixels on a four-unit x/y grid anchored at the",
            "minimum source coordinates; retain original column order."),
        parameters = list(stride = 4L, origin = origin),
        barcodes = colnames(result), features = c(metabolomics = nrow(object)),
        retained = paste("All source-panel features, processed intensities, coordinates, archived SSC",
            "labels and explicitly historical lipid candidates"),
        omitted = "Other historical annotations and derived metadata",
        interpretation = "A spatial subsample of one section; the source panel is already selected")
    S4Vectors::metadata(result)$intensity_provenance <- list(
        assay = "counts",
        processing = "TIC normalization before peak integration in the source workflow",
        full_spectrum_totals = "Unavailable in the bundled object",
        source = paste0("https://github.com/Vitek-Lab/Cardinal3-vignettes/",
            "blob/main/Vignettes/Vignette1/case-study1-vignette.rmd"))
    S4Vectors::metadata(result)$archived_compartments <- list(
        labels = c(`2` = "Adventitia", `5` = "Muscle", `6` = "Urothelium"),
        source = paste0("https://genomicsmachinelearning.github.io/SpaMTP/",
            "articles/Mouse_Urinary_Bladder.html"),
        interpretation = paste("Archived SpaMTP interpretation of SSC clusters;",
            "independent histological validation is required"))
    stopifnot(identical(SummarizedExperiment::assay(result, "counts"),
        SummarizedExperiment::assay(object, "counts")[, keep, drop = FALSE]),
        identical(SpatialExperiment::spatialCoords(result), xy[keep, , drop = FALSE]))
    dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, outputFile, compress = "xz")
    invisible(outputFile)
}

# Retain every pixel from the public pig fetal section and three prespecified
# profile-ion windows used by Cardinal's segmentation tutorial. This helper
# defines a visualization panel, not a peak-picked spectral universe.
preparePigVignetteData <- function(outputFile, nativeDir = NULL,
                                   offline = !is.null(nativeDir)) {
    if (length(outputFile) != 1L || is.na(outputFile) || !nzchar(outputFile) ||
        file.exists(outputFile)) stop("Supply a new output file.")
    object <- SpaMTPData::spaMTPData("pig_fetus_cardinal", version = "1.1.0",
        local_dir = nativeDir, offline = offline)
    stopifnot(identical(dim(object), c(10200L, 4959L)))
    profile <- Cardinal::spectra(object)
    mz <- Cardinal::mz(object)
    xy <- as.matrix(Cardinal::coord(object))[, c("x", "y"), drop = FALSE]
    targets <- c(187.36, 537.08, 840.43)
    tolerance <- 0.5
    windows <- lapply(targets, function(target) which(abs(mz - target) <= tolerance))
    stopifnot(all(lengths(windows) > 0L))
    intensity <- do.call(rbind, lapply(windows, function(i)
        colSums(profile[i, , drop = FALSE])))
    total_ion_current <- colSums(profile)
    stopifnot(all(is.finite(total_ion_current)), all(total_ion_current > 0))
    rownames(intensity) <- sprintf("window_%.2f", targets)
    colnames(intensity) <- paste(xy[, "x"], xy[, "y"], sep = "_")
    rownames(xy) <- colnames(intensity)
    result <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = intensity,
            logcounts = log1p(sweep(intensity, 2, total_ion_current, "/") * 1e4)),
        rowData = S4Vectors::DataFrame(mz = targets, target_mz = targets,
            window_lower = targets - tolerance, window_upper = targets + tolerance,
            source_mz_min = vapply(windows, function(i) min(mz[i]), numeric(1)),
            source_mz_max = vapply(windows, function(i) max(mz[i]), numeric(1)),
            source_bins = lengths(windows),
            source_context = c("Heart-enriched in Cardinal tutorial",
                "Liver-enriched in Cardinal tutorial",
                "Brain/spinal-cord-enriched in Cardinal tutorial"),
            row.names = rownames(intensity)),
        colData = S4Vectors::DataFrame(sample_id = rep("PIGII_206", ncol(object)),
            total_ion_current = total_ion_current, row.names = colnames(intensity)),
        spatialCoords = xy)
    # Cardinal objects do not receive the native-SPE provenance attachment.
    provenance <- list(resource = "pig_fetus_cardinal", version = "1.1.0",
        source_url = "https://zenodo.org/records/17246555",
        md5 = "d056b5c31c136b6e64ffb975ebca9100",
        organism = "Sus scrofa", taxonomy_id = 9823L)
    S4Vectors::metadata(result)$SpaMTPData <- provenance
    S4Vectors::metadata(result)$vignette_subset <- list(
        source = provenance,
        source_url = "https://zenodo.org/records/17246555",
        source_doi = "10.5281/zenodo.17246555",
        source_md5 = "d056b5c31c136b6e64ffb975ebca9100", source_bytes = 9911076L,
        source_license = "CC-BY-4.0", species = "Sus scrofa",
        method = "DESI on a Thermo LTQ instrument",
        selection_source = paste0("https://bioconductor.org/packages/release/data/",
            "experiment/vignettes/CardinalWorkflows/inst/doc/MSI-segmentation.html"),
        original_features = nrow(object), original_observations = ncol(object),
        rule = paste("Retain all pixels. Sum original profile-bin intensities",
            "within plus/minus 0.5 m/z of the three tutorial-specified ions."),
        parameters = list(targets = targets, tolerance_mz = tolerance),
        source_mass_bins = lapply(windows, function(i) mz[i]),
        barcodes = colnames(result),
        normalization = paste("log1p(window intensity / total ion current * 10000);",
            "total ion current uses all 10200 original mass bins"),
        retained = "All 4959 original spatial positions and three profile-ion windows",
        omitted = "All other mass bins; no anatomical labels exist in the source object",
        interpretation = paste("A single fetal section. Organ enrichment is the",
            "original tutorial's interpretation, not pixelwise anatomical truth.",
            "Mass windows are not chemical identifications or a full feature universe."))
    stopifnot(identical(SummarizedExperiment::assay(result, "counts"), intensity),
        identical(SpatialExperiment::spatialCoords(result), xy),
        identical(result$total_ion_current, total_ion_current))
    dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, outputFile, compress = "xz")
    invisible(outputFile)
}
