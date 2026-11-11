# Development assistance: see inst/COPYRIGHTS.

# Small, deterministic teaching data. The identifiers and pathways are fictional.
# Run from the package root after loading SpaMTP.
library(SpaMTP)
metabolites <- data.frame(
    ramp_id = paste0("RAMP_C_", 1:4), chem_source_id = paste0("example:", 1:4),
    chem_data_source = "example",
    common_name = c("2-Hydroxyglutarate", "Glucose", "Lactate", "3-Hydroxyglutarate"),
    mol_formula = c("C5H8O5", "C6H12O6", "C3H6O3", "C5H8O5"),
    iso_smiles = c("OC(CCC(O)=O)C(O)=O", "OCC(O)C(O)C(O)C(O)C=O", "CC(O)C(O)=O",
        "O=C(O)CC(O)CC(O)=O"))
metabolites$monoisotop_mass <- MetaboCoreUtils::calculateMass(metabolites$mol_formula)
metabolites <- annotateSMILESStructure(metabolites)
masses <- c(as.vector(t(MetaboCoreUtils::mass2mz(metabolites$monoisotop_mass[1:3],
    c("[M+H]+", "[M+Na]+")))), 500, 600)
counts <- outer(seq_len(8), seq_len(24), function(i, j)
    as.numeric((i * j + i^2 + j^2) %% 23 + 1))
dimnames(counts) <- list(paste0("ion", seq_len(8)), paste0("pixel", seq_len(24)))
counts[1:2, 1:12] <- counts[1:2, 1:12] + 20
experiment <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts, logcounts = log1p(counts) / log(2)),
    rowData = S4Vectors::DataFrame(mz = masses),
    colData = S4Vectors::DataFrame(region = rep(c("left", "right"), each = 12),
        replicate = rep(paste0("sample", 1:6), each = 4)),
    spatialCoords = cbind(x = rep(1:6, 4), y = rep(1:4, each = 6)))
rna <- counts[1:8, ] + counts[8:1, ]
rownames(rna) <- paste0("GENE", 1:8)
experiment <- addTranscriptome(experiment, rna)
SingleCellExperiment::altExp(experiment, "transcriptome") <- normalizeSMData(
    SingleCellExperiment::altExp(experiment, "transcriptome"),
    "LogNormalize", verbose = FALSE)
index <- buildMZAnnotationIndex(metabolites, adducts = c("M+H", "M+Na"))
experiment <- annotateSM(experiment, index = index, ppm = 5,
    minScore = 0, returnOnlyAnnotated = FALSE, saveIntermediate = TRUE,
    verbose = FALSE)
experiment <- withr::with_seed(1,
    runMetabolicPCA(experiment, npcs = 3, verbose = FALSE))
experiment <- experiment[order(SummarizedExperiment::rowData(experiment)$mz), ]
gene_reference <- data.frame(hgnc_id = paste0("HGNC:", 1:8),
    symbol = paste0("GENE", 1:8), status = "Approved",
    entrez_id = as.character(101:108), ensembl_gene_id = paste0("ENSG000", 1:8),
    uniprot_ids = paste0("U", 1:8), prev_symbol = c("OLD1", rep("", 7)))
source <- data.frame(rampId = c(paste0("RAMP_G_", 1:8), metabolites$ramp_id),
    sourceId = c(paste0("entrez:", 101:108), metabolites$chem_source_id),
    commonName = c(paste0("GENE", 1:8), metabolites$common_name))
sets <- list(P1 = c(paste0("RAMP_G_", 1:4), "RAMP_C_1", "RAMP_C_2"),
    P2 = c(paste0("RAMP_G_", 4:7), "RAMP_C_2", "RAMP_C_3"),
    P3 = c(paste0("RAMP_G_", c(2, 5, 8)), "RAMP_C_3", "RAMP_C_4"))
database <- list(chem_props = metabolites, source_df = source,
    gene_reference = gene_reference,
    analyte = data.frame(rampId = source$rampId),
    analytehaspathway = data.frame(rampId = unlist(sets, use.names = FALSE),
        pathwayRampId = rep(names(sets), lengths(sets))),
    pathway = data.frame(pathwayRampId = names(sets), sourceId = names(sets),
        type = "example",
        pathwayName = c("Example pathway 1", "Example pathway 2", "Example pathway 3")))
spamtp_example <- list(experiment = experiment, metabolites = metabolites,
    database = database, gene_reference = gene_reference)
save(spamtp_example, file = "data/spamtp_example.rda", compress = "xz", version = 2)
