.dehydration_test_mass <- function(carbon, hydrogen, oxygen, nitrogen = 0) {
  12 * carbon + 1.00782503223 * hydrogen +
    15.99491461957 * oxygen + 14.00307400443 * nitrogen
}

.dehydration_test_sphingosine <- function() {
  data.frame(
    name = "sphingosine", formula = "C18H37NO2",
    exactmass = .dehydration_test_mass(18, 37, 2, 1),
    iso_smiles = "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](N)CO"
  )
}

.dehydration_test_mz <- function(neutral_mass, losses) {
  neutral_mass + 1.007276466621 - losses * (2 * 1.00782503223 + 15.99491461957)
}

test_that("water-loss rules preserve charge and use neutral water masses", {
  rules <- adductRules("positive")
  selected <- match(c("M+H-H2O", "M+H-2H2O"), rules$name)
  expect_false(anyNA(selected))
  losses <- rules[selected, , drop = FALSE]
  expect_equal(losses$mass_shift, .dehydration_test_mz(0, 1:2),
               tolerance = 1e-10)
  expect_equal(losses$n_molecules, c(1, 1))
  expect_equal(losses$charge, c(1, 1))
  expect_equal(losses$loss_h, c(0, 0))
  expect_equal(losses$water_loss, 1:2)
  expect_equal(losses$added_charge - losses$loss_h, losses$charge)
  expect_true(all(losses$prior < rules$prior[rules$name == "M+H"]))
  expect_false(any(losses$name %in% adductRules("negative")$name))
  expect_true(all(losses$name %in% maldiMatrixRules("DHB")$name))
})

test_that("sphingosine produces structurally supported single and double water loss", {
  target <- .dehydration_test_sphingosine()
  observed <- .dehydration_test_mz(target$exactmass, 1:2)
  index <- buildMZAnnotationIndex(
    target, polarity = "positive", adducts = c("M+H-H2O", "M+H-2H2O")
  )
  result <- queryMZAnnotationIndex(observed, index, ppm = 1e-6)
  expect_equal(nrow(result), 2)
  expect_setequal(result$adduct, c("M+H-H2O", "M+H-2H2O"))
  expect_true(all(result$charge == 1))
  expect_true(all(result$metabolite_names == "sphingosine"))
  expect_true(all(result$max_water_losses >= result$water_loss))
  expect_equal(sort(result$expected_mz), sort(observed), tolerance = 1e-10)

  measured <- queryMZAnnotationIndex(264.268615, index, ppm = 0.2)
  expect_equal(measured$adduct, "M+H-2H2O")
  expect_lt(measured$ppm_error, 0.2)
  expect_length(queryMZAnnotationIndex(264.268615, index, ppm = 0.01)$adduct, 0)
})

test_that("oxygen presence alone does not admit water-loss candidates", {
  target <- data.frame(
    name = c("dimethyl ether", "phenol", "acetic acid", "neopentyl alcohol",
             "methanol"),
    formula = c("C2H6O", "C6H6O", "C2H4O2", "C5H12O", "CH4O"),
    exactmass = .dehydration_test_mass(
      c(2, 6, 2, 5, 1), c(6, 6, 4, 12, 4), c(1, 1, 2, 1, 1)),
    iso_smiles = c("COC", "Oc1ccccc1", "CC(=O)O", "CC(C)(C)CO", "CO")
  )
  for (i in seq_len(nrow(target))) {
    index <- buildMZAnnotationIndex(
      target[i, , drop = FALSE], polarity = "positive",
      adducts = c("M+H-H2O", "M+H-2H2O"), min_structure_score = 0
    )
    observed <- .dehydration_test_mz(target$exactmass[[i]], 1:2)
    observed <- observed[observed > 0]
    expect_equal(nrow(queryMZAnnotationIndex(observed, index, ppm = 1e-6)), 0,
                 info = target$name[[i]])
  }
})

test_that("unknown or invalid structures cannot bypass water-loss eligibility", {
  target <- data.frame(
    name = c("missing", "empty", "invalid"), formula = "C2H6O",
    exactmass = .dehydration_test_mass(2, 6, 1),
    iso_smiles = c(NA_character_, "", "C1CC")
  )
  index <- buildMZAnnotationIndex(
    target, adducts = c("M+H", "M+H-H2O", "M+H-2H2O"),
    collapse_isomers = FALSE, min_structure_score = 0
  )
  result <- queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass[[1]], 0:2), index, ppm = 1e-6
  )
  expect_setequal(result$metabolite_names, target$name)
  expect_true(all(result$adduct == "M+H"))

  target$iso_smiles <- NULL
  no_smiles <- buildMZAnnotationIndex(
    target, adducts = c("M+H-H2O", "M+H-2H2O"),
    collapse_isomers = FALSE, min_structure_score = 0
  )
  expect_equal(nrow(queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass[[1]], 1:2), no_smiles, ppm = 1e-6)), 0)
})

test_that("two hydroxyl groups cannot reuse a single neighbouring hydrogen", {
  structures <- deconvolveSMILES(c(
    one_shared_hydrogen = "OCC(C)CO",
    two_shared_hydrogens = "OCCCO"
  ), strict = TRUE)
  expect_equal(structures$dehydration_sites, c(2, 2))
  expect_equal(structures$max_water_losses, c(1, 2))

  predictions <- lapply(structures$smiles, function(smiles) {
    predictAdductsFromSMILES(smiles, polarity = "positive")
  })
  expect_true(predictions[[1]]$retained[
    predictions[[1]]$adduct == "M+H-H2O"])
  expect_false(predictions[[1]]$retained[
    predictions[[1]]$adduct == "M+H-2H2O"])
  expect_true(predictions[[2]]$retained[
    predictions[[2]]$adduct == "M+H-2H2O"])
})

test_that("water losses stay within one connected analyte and respect explicit H", {
  structures <- deconvolveSMILES(c(
    separate_ethanols = "CCO.CCO",
    explicit_hydrogens = "[H]OC([H])([H])C([H])([H])[H]",
    enol = "C=C(O)C"
  ), strict = TRUE)
  expect_equal(structures$dehydration_sites, c(2, 1, 0))
  expect_equal(structures$max_water_losses, c(1, 1, 0))
})

test_that("explicit adduct subsets continue to bound the search space", {
  target <- .dehydration_test_sphingosine()
  protonated <- buildMZAnnotationIndex(target, adducts = "M+H")
  result <- queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass, 0:2), protonated, ppm = 1e-6
  )
  expect_equal(result$adduct, "M+H")
  expect_equal(protonated$rules$name, "M+H")

  double_loss <- buildMZAnnotationIndex(target, adducts = "[M+H-2H2O]+")
  result <- queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass, 0:2), double_loss, ppm = 1e-6
  )
  expect_equal(result$adduct, "M+H-2H2O")
})

test_that("annotation wrappers preserve gated water-loss candidates", {
  target <- .dehydration_test_sphingosine()
  target$id <- "sphingosine"
  observed <- .dehydration_test_mz(target$exactmass, 2)
  index <- buildMZAnnotationIndex(target, adducts = "M+H-2H2O")
  direct <- annotateMZ(observed, index = index, ppm = 1e-6)
  expect_equal(direct$adduct, "M+H-2H2O")
  batch <- annotateBigData(
    c(observed, 999), index = index, ppm_error = 1e-6, verbose = FALSE
  )
  expect_equal(nrow(batch), 1)
  expect_match(batch$all_Adducts, "M+H-2H2O", fixed = TRUE)

  counts <- matrix(c(10, 3, 20, 0), nrow = 2,
                   dimnames = list(c("dehydrated", "unmatched"), c("p1", "p2")))
  experiment <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(mz = c(observed, 999)),
    spatialCoords = matrix(c(0, 1, 0, 0), nrow = 2,
                           dimnames = list(c("p1", "p2"), c("x", "y")))
  )
  annotated <- annotateSM(
    experiment, index = index, raw.mz.column = "mz", ppm_error = 1e-6,
    return.only.annotated = FALSE, verbose = FALSE
  )
  expect_identical(SummarizedExperiment::assay(annotated, "counts"), counts)
  candidates <- S4Vectors::metadata(annotated)$mz_annotation$results
  expect_equal(candidates$Adduct, "M+H-2H2O")
})

test_that("malformed custom water-loss rules fail before candidate search", {
  target <- .dehydration_test_sphingosine()
  rule <- adductRules("positive")
  rule <- rule[rule$name == "M+H-2H2O", , drop = FALSE]
  bad_mass <- rule
  bad_mass$mass_shift <- bad_mass$mass_shift + 0.01
  expect_error(buildMZAnnotationIndex(target, rules = bad_mass), "Water-loss")
  bad_charge <- rule
  bad_charge$charge <- bad_charge$added_charge <- 2L
  expect_error(buildMZAnnotationIndex(target, rules = bad_charge), "Water-loss")
  for (loss in c(-1, 1.5, 3, NA_real_)) {
    bad_count <- rule
    bad_count$water_loss <- loss
    expect_error(buildMZAnnotationIndex(target, rules = bad_count), "water_loss")
  }
  bad_name <- rule
  bad_name$water_loss <- 0L
  expect_error(buildMZAnnotationIndex(target, rules = bad_name),
               "[Ww]ater.loss|water_loss")
  bad_name$name <- "[M+H-2H2O]+"
  expect_error(buildMZAnnotationIndex(target, rules = bad_name),
               "[Ww]ater.loss|water_loss")

  old_schema <- rule
  old_schema$water_loss <- NULL
  target$iso_smiles <- NULL
  index <- buildMZAnnotationIndex(target, rules = old_schema)
  expect_equal(index$rules$water_loss, 2L)
  expect_equal(nrow(queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass, 2), index, ppm = 1e-6)), 0)
})

test_that("legacy structural tables parse only mass-matched water-loss candidates", {
  target <- rbind(
    .dehydration_test_sphingosine(),
    data.frame(
      name = c("ethanol", "phenol", "acetic acid"),
      formula = c("C2H6O", "C6H6O", "C2H4O2"),
      exactmass = .dehydration_test_mass(c(2, 6, 2), c(6, 6, 4), c(1, 1, 2)),
      iso_smiles = c("CCO", "Oc1ccccc1", "CC(=O)O")
    )
  )
  legacy <- annotateSMILESStructure(target)
  legacy$dehydration_sites <- NULL
  legacy$max_water_losses <- NULL
  withr::local_options(SpaMTP.max_runtime_smiles = 0L)

  parsed_smiles <- character()
  original_parser <- getFromNamespace(".parse_smiles_graph", "SpaMTP")
  local_mocked_bindings(
    .parse_smiles_graph = function(smiles) {
      parsed_smiles <<- c(parsed_smiles, smiles)
      original_parser(smiles)
    },
    .package = "SpaMTP"
  )
  index <- buildMZAnnotationIndex(
    legacy, adducts = c("M+H", "M+H-H2O", "M+H-2H2O"),
    collapse_isomers = FALSE
  )
  expect_length(parsed_smiles, 0)
  result <- queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass[[1]], 2), index, ppm = 1e-6
  )
  expect_equal(result$metabolite_names, "sphingosine")
  expect_equal(result$adduct, "M+H-2H2O")
  expect_setequal(parsed_smiles, target$iso_smiles[[1]])
  expect_length(parsed_smiles, 1)
  unsupported <- .dehydration_test_mz(target$exactmass[3:4], 1)
  expect_equal(nrow(queryMZAnnotationIndex(unsupported, index, ppm = 1e-6)), 0)
  expect_length(parsed_smiles, 1)
})

test_that("legacy isomer collapsing preserves distinct water-loss eligibility", {
  target <- data.frame(
    name = c("shared hydrogen", "independent hydrogens"),
    formula = "C4H10O2",
    exactmass = .dehydration_test_mass(4, 10, 2),
    iso_smiles = c("OCC(C)CO", "OCCC(O)C")
  )
  legacy <- annotateSMILESStructure(target)
  legacy$dehydration_sites <- NULL
  legacy$max_water_losses <- NULL
  index <- buildMZAnnotationIndex(
    legacy, adducts = "M+H-2H2O", collapse_isomers = TRUE
  )
  result <- queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass[[1]], 2), index, ppm = 1e-6
  )
  expect_equal(result$metabolite_names, "independent hydrogens")
})

test_that("infer_structure never uses only existing water-loss eligibility", {
  target <- .dehydration_test_sphingosine()
  precomputed <- annotateSMILESStructure(target)
  local_mocked_bindings(
    .smiles_dehydration_capacity = function(...) {
      stop("infer_structure = 'never' must not invoke the capacity parser")
    },
    .package = "SpaMTP"
  )
  unknown <- buildMZAnnotationIndex(
    target, adducts = "M+H-2H2O", infer_structure = "never"
  )
  known <- buildMZAnnotationIndex(
    precomputed, adducts = "M+H-2H2O", infer_structure = "never"
  )
  observed <- .dehydration_test_mz(target$exactmass, 2)
  expect_equal(nrow(queryMZAnnotationIndex(observed, unknown, ppm = 1e-6)), 0)
  expect_equal(
    queryMZAnnotationIndex(observed, known, ppm = 1e-6)$metabolite_names,
    "sphingosine"
  )
})

test_that("lazy structure resolution preserves complete scoring and score cutoffs", {
  target <- .dehydration_test_sphingosine()
  precomputed <- annotateSMILESStructure(target)
  complete <- buildMZAnnotationIndex(
    precomputed, adducts = "M+H-2H2O", min_structure_score = 0.5
  )
  withr::local_options(SpaMTP.max_runtime_smiles = 0L)
  expect_warning(
    lazy <- buildMZAnnotationIndex(
      target, adducts = "M+H-2H2O", min_structure_score = 0.5
    ),
    "Skipping runtime SMILES decomposition"
  )
  observed <- .dehydration_test_mz(target$exactmass, 2)
  expected <- queryMZAnnotationIndex(observed, complete, ppm = 1e-6)
  actual <- queryMZAnnotationIndex(observed, lazy, ppm = 1e-6)
  expect_equal(nrow(expected), 1)
  expect_equal(nrow(actual), 1)
  columns <- c("adduct", "structure_score", "structure_status", "score",
               "dehydration_sites", "max_water_losses", "dehydration_status")
  expect_equal(actual[, columns], expected[, columns], tolerance = 1e-12)
})

test_that("adding water-loss rules preserves results for existing adducts", {
  target <- .dehydration_test_sphingosine()
  all_rules <- adductRules("positive")
  old_rules <- all_rules[
    !all_rules$name %in% c("M+H-H2O", "M+H-2H2O"), , drop = FALSE]
  before <- buildMZAnnotationIndex(target, rules = old_rules)
  after <- buildMZAnnotationIndex(target)
  observed <- before$expected_mz
  baseline <- queryMZAnnotationIndex(observed, before, ppm = 1e-6)
  current <- queryMZAnnotationIndex(observed, after, ppm = 1e-6)
  current <- current[!current$adduct %in% c("M+H-H2O", "M+H-2H2O"), ]
  columns <- c("observed_mz", "expected_mz", "ppm_error", "adduct", "charge",
               "metabolite_names", "mass_score", "rule_prior", "chemical_score",
               "structure_score", "reactive_site_score", "score")
  rownames(current) <- NULL
  expect_equal(current[, columns], baseline[, columns], tolerance = 1e-12)
})

test_that("ordinary adducts retain duplicate records when collapsing is disabled", {
  target <- .dehydration_test_sphingosine()
  target <- rbind(target, target)
  index <- buildMZAnnotationIndex(
    target, adducts = c("M+H", "M+H-2H2O"), collapse_isomers = FALSE
  )
  result <- queryMZAnnotationIndex(
    .dehydration_test_mz(target$exactmass[[1]], 0), index, ppm = 1e-6
  )
  expect_equal(nrow(result), 2)
  expect_true(all(result$adduct == "M+H"))
})
