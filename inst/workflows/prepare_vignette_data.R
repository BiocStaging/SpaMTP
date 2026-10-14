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
            retained = "All measured features, raw counts, pairing and original coordinates",
            omitted = "Archived annotations, normalized layers, reductions and images",
            interpretation = "A spatial subset of one specimen; not biological replication")
        result
    }
    dhb <- read("mouse_brain_dhb_striatum")
    xy <- SpatialExperiment::spatialCoords(dhb)
    keep <- lapply(c("intact", "lesioned"), function(side) {
        center <- colMeans(xy[dhb$RegionLoupe == "ACB" & dhb$lesion == side, ,
                              drop = FALSE])
        eligible <- which(dhb$lesion == side)
        distance <- rowSums(sweep(xy[eligible, , drop = FALSE], 2, center)^2)
        eligible[order(distance, colnames(dhb)[eligible])[seq_len(80L)]]
    })
    result <- subset(dhb, sort(unique(unlist(keep))),
        c("sample_id", "RegionLoupe", "lesion"),
        paste("Within each lesion side, take the 80 locations nearest the",
              "centroid of archived ACB coordinates; break distance ties by barcode",
              "and retain the original column order."),
        list(per_side = 80L, center_region = "ACB", sides = c("intact", "lesioned")))
    saveRDS(result, file.path(outputDir, "dhb_demo.rds"), compress = "xz")
    human <- read("human_brain_integrated")
    xy <- SpatialExperiment::spatialCoords(human)
    center <- apply(xy, 2, stats::median)
    distance <- rowSums(sweep(xy, 2, center)^2)
    keep <- sort(order(distance, colnames(human))[seq_len(240L)])
    result <- subset(human, keep, c("sample_id", "cout_labels"),
        paste("Take the 240 locations nearest the coordinate-wise median;",
              "break distance ties by barcode and retain original column order."),
        list(observations = 240L, center = center))
    saveRDS(result, file.path(outputDir, "human_demo.rds"), compress = "xz")
    invisible(outputDir)
}
