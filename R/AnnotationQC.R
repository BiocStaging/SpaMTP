## Quality control for MS1 m/z annotation: mass-shift null diagnostics,
## isotope and matrix peaks, MALDI-plausible rule sets and adduct groups.
## All functions are opt-in; they do not change annotateSM() defaults.

#' Assess annotation matches with mass-shift decoys
#'
#' Scores the observed peak list and decoy copies of it in which every peak,
#' and the optional MS1 context spectrum, is moved by the same mass (in Da).
#' A constant shift larger than twice the tolerance moves every peak away
#' from its own true identity while keeping isotope and adduct spacings
#' exact. The decoy-to-target ratio of annotated peaks at each score cutoff
#' measures how often perturbed spectra still obtain matches. It is an
#' exploratory specificity diagnostic, not a calibrated probability that a
#' reported metabolite identity is wrong.
#'
#' Interpreting the ratio as an FDR estimate requires decoy matches to model
#' incorrect target matches. That assumption is not established by shifting
#' masses alone: shift magnitude can change local database density, mass
#' defects and contextual adduct evidence. Compare several positive and
#' negative shifts and retain the target counts and per-shift scores. A low
#' ratio, especially with few retained peaks or few decoys, does not by itself
#' demonstrate FDR control. Chemical identity requires independent evidence.
#'
#' @param observed_mz Numeric vector of observed m/z values (one per peak).
#' @param index A `spamtp_mz_index` from [buildMZAnnotationIndex()].
#' @param ppm Mass tolerance in parts per million.
#' @param ms1_spectrum Optional data frame with `mz` and `intensity` columns
#'   used as isotope and adduct-family context; it is shifted with the peaks.
#' @param shifts_da Decoy shifts in Da. Each must exceed twice the tolerance
#'   at the largest observed m/z. By default, -8, -5, 5 and 8 times that
#'   tolerance half-width (`ppm` at the largest m/z).
#' @param cutoffs Score cutoffs at which to report the estimate.
#' @param ... Further arguments passed to [queryMZAnnotationIndex()].
#'
#' @return An object of class `spamtp_annotation_fdr`: a list with `curve`
#'   (cutoff, annotated target peaks, mean annotated decoy peaks and the
#'   decoy-to-target ratio, capped at 1, in the historical `fdr` column),
#'   `best_scores` (best score per peak for the target and each decoy;
#'   0 means no candidate within tolerance) and the settings.
#' @seealso [annotationScoreCutoff()]
#' @examples
#' db <- data.frame(
#'   formula = c("C5H8O5", "C6H12O6", "C4H9NO2"),
#'   exactmass = c(148.037173366, 180.063388104, 103.063328537),
#'   isomers_names = c("2-hydroxyglutarate", "hexose", "GABA")
#' )
#' index <- buildMZAnnotationIndex(db, adducts = c("M+H", "M+Na", "M+K"))
#' fdr <- estimateAnnotationFDR(c(149.04445, 203.05261, 250.1), index)
#' fdr$curve[fdr$curve$cutoff %in% c(0, 0.5), ]
#' @export
estimateAnnotationFDR <- function(observed_mz, index, ppm = 5,
                                  ms1_spectrum = NULL, shifts_da = NULL,
                                  cutoffs = seq(0, 1, by = 0.01), ...) {
  if (!inherits(index, "spamtp_mz_index")) {
    stop("index must be created by buildMZAnnotationIndex().", call. = FALSE)
  }
  observed_mz <- suppressWarnings(as.numeric(observed_mz))
  if (!length(observed_mz) || any(!is.finite(observed_mz) | observed_mz <= 0)) {
    stop("observed_mz must contain positive finite values.", call. = FALSE)
  }
  half_width <- ppm * 1e-6 * max(observed_mz)
  if (is.null(shifts_da)) shifts_da <- c(-8, -5, 5, 8) * half_width
  if (!is.numeric(shifts_da) || !length(shifts_da) ||
      any(!is.finite(shifts_da)) || any(abs(shifts_da) <= 2 * half_width)) {
    stop("Each decoy shift must exceed twice the tolerance at the largest ",
         "m/z (", signif(2 * half_width, 3), " Da).", call. = FALSE)
  }
  if (!is.null(ms1_spectrum) &&
      !all(c("mz", "intensity") %in% names(ms1_spectrum))) {
    stop("ms1_spectrum must have mz and intensity columns.", call. = FALSE)
  }

  best_for <- function(shift) {
    query <- observed_mz + shift
    spectrum <- if (is.null(ms1_spectrum)) NULL else data.frame(
      mz = ms1_spectrum$mz + shift, intensity = ms1_spectrum$intensity
    )
    res <- queryMZAnnotationIndex(query, index, ppm = ppm,
                                  ms1_spectrum = spectrum, min_score = 0, ...)
    values <- unique(query)
    best <- numeric(length(values))
    if (nrow(res)) {
      agg <- tapply(res$score, match(res$observed_mz, values), max)
      best[as.integer(names(agg))] <- agg
    }
    best[match(query, values)]
  }

  sets <- c(0, shifts_da)
  best <- vapply(sets, best_for, numeric(length(observed_mz)))
  best <- matrix(best, nrow = length(observed_mz))
  colnames(best) <- c("target", sprintf("decoy_%+.4gDa", shifts_da))
  annotated <- function(scores, cutoff) {
    if (cutoff <= 0) sum(scores > 0) else sum(scores >= cutoff)
  }
  curve <- do.call(rbind, lapply(cutoffs, function(cutoff) {
    target <- annotated(best[, 1], cutoff)
    decoy <- mean(apply(best[, -1, drop = FALSE], 2, annotated, cutoff))
    data.frame(cutoff = cutoff, target_peaks = target, decoy_peaks = decoy,
               fdr = if (target > 0) min(1, decoy / target) else NA_real_)
  }))
  structure(
    list(curve = curve,
         best_scores = data.frame(mz = observed_mz, best, check.names = FALSE),
         ppm = ppm, shifts_da = shifts_da),
    class = "spamtp_annotation_fdr"
  )
}

#' @rdname estimateAnnotationFDR
#' @param x A `spamtp_annotation_fdr` object.
#' @method print spamtp_annotation_fdr
#' @export
print.spamtp_annotation_fdr <- function(x, ...) {
  curve <- x$curve
  show <- curve[curve$cutoff %in% intersect(c(0, 0.05, 0.1, 0.2, 0.5),
                                            curve$cutoff), , drop = FALSE]
  cat("SpaMTP annotation FDR (mass-shift decoys:",
      paste0(sprintf("%+.4g", x$shifts_da), collapse = ", "), "Da; tolerance",
      x$ppm, "ppm)\n")
  print(show, row.names = FALSE)
  invisible(x)
}

#' Choose a cutoff from the mass-shift diagnostic
#'
#' Finds a cutoff meeting a requested decoy-to-target match ratio. The
#' result is a diagnostic threshold; it does not guarantee control of the
#' false-discovery rate. See the assumptions in [estimateAnnotationFDR()].
#'
#' @param fdr A `spamtp_annotation_fdr` object from [estimateAnnotationFDR()].
#' @param target_fdr Largest acceptable decoy-to-target match ratio.
#' @param min_peaks Smallest number of annotated target peaks a cutoff must
#'   keep; estimates from a handful of peaks are too noisy to act on.
#'
#' @return The smallest score cutoff whose decoy-to-target ratio does not exceed
#'   `target_fdr` while keeping at least `min_peaks` annotated peaks, or `NA`
#'   when no cutoff qualifies.
#' @examples
#' fdr <- structure(list(curve = data.frame(
#'   cutoff = c(0, 0.1, 0.5), target_peaks = c(100, 60, 10),
#'   decoy_peaks = c(90, 20, 1), fdr = c(0.9, 1 / 3, 0.1))),
#'   class = "spamtp_annotation_fdr")
#' annotationScoreCutoff(fdr, target_fdr = 0.2)
#' @export
annotationScoreCutoff <- function(fdr, target_fdr = 0.1, min_peaks = 10) {
  if (!inherits(fdr, "spamtp_annotation_fdr")) {
    stop("fdr must come from estimateAnnotationFDR().", call. = FALSE)
  }
  ok <- is.finite(fdr$curve$fdr) & fdr$curve$fdr <= target_fdr &
    fdr$curve$target_peaks >= max(1, min_peaks)
  if (!any(ok)) return(NA_real_)
  min(fdr$curve$cutoff[ok])
}

#' MALDI-plausible adduct rules
#'
#' Returns a restricted matrix rule set for the `adducts` argument of
#' [buildMZAnnotationIndex()]. By default it retains singly charged rules and
#' removes the solvent families listed in the implementation. Select rules
#' using the acquisition protocol and experimental evidence; a preset does
#' not establish that excluded ions cannot occur or that retained matches
#' are correct.
#'
#' @param maldi_matrix MALDI matrix name, as in [maldiMatrixRules()].
#' @param polarity Optional polarity; defaults to the matrix's polarity.
#' @param singly_charged Keep only singly charged rules.
#' @param exclude_solvent Drop electrospray solvent and solvent-multimer rules.
#' @param include_matrix_adducts Keep matrix-specific adduct rules.
#'
#' @return A character vector of rule names.
#' @examples
#' maldiAdductPreset("DHB")
#' @export
maldiAdductPreset <- function(maldi_matrix, polarity = NULL,
                              singly_charged = TRUE, exclude_solvent = TRUE,
                              include_matrix_adducts = TRUE) {
  rules <- maldiMatrixRules(maldi_matrix, polarity = polarity)
  keep <- rep(TRUE, nrow(rules))
  if (isTRUE(singly_charged)) keep <- keep & abs(rules$charge) == 1
  if (isTRUE(exclude_solvent)) {
    solvent <- grepl("ACN|CH3OH|IsoProp|DMSO|FA-|Hac|TFA", rules$name)
    if ("complexity" %in% names(rules)) {
      solvent <- solvent | rules$complexity %in% "solvent"
    }
    keep <- keep & !solvent
  }
  if (!isTRUE(include_matrix_adducts) && "rule_class" %in% names(rules)) {
    keep <- keep & !rules$rule_class %in% "matrix_adduct"
  }
  unique(rules$name[keep])
}

.matrix_formulas <- function() {
  c(dhb = "C7H6O4", chca = "C10H7NO3", `9aa` = "C13H10N2",
    dan = "C10H10N2", norharmane = "C11H8N2")
}

.formula_mass <- function(formula) {
  m <- .spamtp_mass_constants()
  mass <- c(C = m$carbon, H = m$hydrogen, N = m$nitrogen, O = m$oxygen)
  sum(vapply(names(mass), function(el) {
    .element_count(formula, el) * mass[[el]]
  }, numeric(1)))
}

#' Ions formed by the MALDI matrix itself
#'
#' Enumerates matrix monomer and cluster ions: protonated, sodiated and
#' potassiated n-mers with up to one water loss per molecule, metal-exchanged
#' n-mers and the radical cation in positive mode, and deprotonated n-mers
#' with water losses in negative mode.
#'
#' @param maldi_matrix One of `"dhb"`, `"chca"`, `"9aa"`, `"dan"` or
#'   `"norharmane"` (case-insensitive aliases accepted).
#' @param polarity `"positive"` or `"negative"`.
#' @param max_molecules Largest cluster size.
#'
#' @return A data frame with the ion label, cluster size, water losses and m/z.
#' @examples
#' head(matrixIonTable("DHB"))
#' @export
matrixIonTable <- function(maldi_matrix = "dhb",
                           polarity = c("positive", "negative"),
                           max_molecules = 3L) {
  polarity <- match.arg(polarity)
  key <- .normalise_maldi_matrix(maldi_matrix)
  formulas <- .matrix_formulas()
  if (!key %in% names(formulas)) {
    stop("No matrix ion formula is registered for '", maldi_matrix, "'.",
         call. = FALSE)
  }
  m <- .spamtp_mass_constants()
  formula <- formulas[[key]]
  mass <- .formula_mass(formula)
  has_o <- .element_count(formula, "O") > 0
  label <- toupper(key)
  rows <- list()
  add <- function(n, k, name, shift) {
    rows[[length(rows) + 1L]] <<- data.frame(
      matrix = label, ion = name, n_molecules = n, water_loss = k,
      mz = n * mass - k * m$water + shift, stringsAsFactors = FALSE)
  }
  for (n in seq_len(max_molecules)) {
    mm <- if (n == 1) "M" else paste0(n, "M")
    losses <- if (has_o) 0:n else 0L
    for (k in losses) {
      w <- if (k == 0) "" else if (k == 1) "-H2O" else paste0("-", k, "H2O")
      if (polarity == "positive") {
        add(n, k, paste0("[", mm, "+H", w, "]+"), m$proton)
        add(n, k, paste0("[", mm, "+Na", w, "]+"), m$sodium_ion)
        add(n, k, paste0("[", mm, "+K", w, "]+"), m$potassium_ion)
      } else {
        add(n, k, paste0("[", mm, "-H", w, "]-"), -m$proton)
      }
    }
    if (polarity == "positive") {
      add(n, 0L, paste0("[", mm, "+2Na-H]+"), 2 * m$sodium_ion - m$proton)
      add(n, 0L, paste0("[", mm, "+2K-H]+"), 2 * m$potassium_ion - m$proton)
    }
  }
  if (polarity == "positive") add(1L, 0L, "[M]+.", -m$electron)
  out <- do.call(rbind, rows)
  out[order(out$mz), , drop = FALSE]
}

#' Flag peaks explained by MALDI matrix ions
#'
#' @param mz Numeric vector of observed m/z values.
#' @param maldi_matrix Matrix name, as in [matrixIonTable()].
#' @param polarity `"positive"` or `"negative"`.
#' @param ppm Mass tolerance in parts per million.
#' @param max_molecules Largest matrix cluster considered.
#'
#' @return A data frame with one row per peak: `matrix_ion` (logical), the
#'   matched ion label and the signed ppm error.
#' @examples
#' flagMatrixPeaks(c(137.0238, 155.0344, 200.1), "DHB")
#' @export
flagMatrixPeaks <- function(mz, maldi_matrix = "dhb",
                            polarity = c("positive", "negative"), ppm = 5,
                            max_molecules = 3L) {
  mz <- suppressWarnings(as.numeric(mz))
  if (any(!is.finite(mz))) stop("mz must be finite.", call. = FALSE)
  ions <- matrixIonTable(maldi_matrix, match.arg(polarity), max_molecules)
  hit <- vapply(mz, function(x) {
    err <- (x - ions$mz) / ions$mz * 1e6
    j <- which.min(abs(err))
    if (abs(err[j]) <= ppm) j else NA_integer_
  }, integer(1))
  data.frame(mz = mz, matrix_ion = !is.na(hit),
             ion = ifelse(is.na(hit), NA_character_, ions$ion[hit]),
             ppm_error = ifelse(is.na(hit), NA_real_,
                                (mz - ions$mz[hit]) / ions$mz[hit] * 1e6),
             stringsAsFactors = FALSE)
}

.feature_profiles <- function(mz, intensity) {
  if (is.null(dim(intensity))) {
    intensity <- as.numeric(intensity)
    if (length(intensity) != length(mz)) {
      stop("intensity must have one value per m/z.", call. = FALSE)
    }
    profiles <- NULL
  } else {
    if (nrow(intensity) != length(mz)) {
      stop("intensity must have one row per m/z (features x observations).",
           call. = FALSE)
    }
    profiles <- as.matrix(intensity)
    intensity <- rowMeans(profiles)
  }
  if (any(!is.finite(intensity)) ||
      (!is.null(profiles) && any(!is.finite(profiles)))) {
    stop("intensity must be finite.", call. = FALSE)
  }
  if (any(intensity < 0)) {
    stop("mean intensities must be non-negative.", call. = FALSE)
  }
  list(mean = intensity, profiles = profiles)
}

.exclusion_mask <- function(exclude, n) {
  if (is.null(exclude)) return(rep(FALSE, n))
  if (is.logical(exclude)) {
    if (length(exclude) != n) {
      stop("A logical exclude must have one value per m/z.", call. = FALSE)
    }
    return(exclude %in% TRUE)
  }
  idx <- suppressWarnings(as.integer(exclude))
  if (any(is.na(idx) | idx < 1L | idx > n)) {
    stop("exclude must be logical or valid peak indices.", call. = FALSE)
  }
  seq_len(n) %in% idx
}

.profile_cor <- function(profiles, i, j) {
  if (is.null(profiles)) return(NA_real_)
  a <- profiles[i, ]
  b <- profiles[j, ]
  if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  stats::cor(a, b)
}

#' Flag peaks that are 13C isotopologues of other peaks
#'
#' A peak is flagged when another peak sits one 13C spacing below it within
#' `ppm`, the heavier peak is weaker than the carbon count of the lighter ion
#' allows, and, when spatial profiles are given, the two co-localise. Flagged
#' peaks should not be annotated as monoisotopic ions.
#'
#' @param mz Numeric vector of observed m/z values.
#' @param intensity Either one mean intensity per peak, or a features x
#'   observations matrix (rows in the order of `mz`); a matrix also enables
#'   the co-localisation test.
#' @param ppm Spacing tolerance in parts per million of the heavier peak.
#' @param charge Ion charge (spacing is 1.00335 / charge).
#' @param min_cor Smallest profile correlation accepted when a matrix is given.
#' @param max_ratio Optional fixed upper bound on heavy/light intensity. By
#'   default it is 1.2 x 1.07% per carbon, with at most m/z x charge / 12
#'   carbons in the lighter ion.
#'
#' @return A data frame with one row per peak: `isotope` (logical), the
#'   index and m/z of the monoisotopic parent, the spacing error, the
#'   intensity ratio, its bound and the profile correlation.
#' @examples
#' flagIsotopePeaks(c(399.1436, 400.1466, 410.2), c(100, 14, 50))
#' @export
flagIsotopePeaks <- function(mz, intensity, ppm = 3, charge = 1L,
                             min_cor = 0.5, max_ratio = NULL) {
  mz <- suppressWarnings(as.numeric(mz))
  if (any(!is.finite(mz) | mz <= 0)) {
    stop("mz must contain positive finite values.", call. = FALSE)
  }
  prof <- .feature_profiles(mz, intensity)
  spacing <- .spamtp_mass_constants()$carbon13_spacing / abs(charge)
  ord <- order(mz)
  smz <- mz[ord]
  out <- data.frame(mz = mz, isotope = FALSE, parent = NA_integer_,
                    parent_mz = NA_real_, ppm_error = NA_real_,
                    intensity_ratio = NA_real_, max_ratio = NA_real_,
                    correlation = NA_real_)
  for (j in seq_along(mz)) {
    target <- mz[j] - spacing
    tol <- mz[j] * ppm * 1e-6
    lo <- findInterval(target - tol, smz) + 1L
    hi <- findInterval(target + tol, smz)
    if (hi < lo) next
    best <- NULL
    for (i in ord[lo:hi]) {
      if (prof$mean[i] <= 0) next
      ratio <- prof$mean[j] / prof$mean[i]
      bound <- if (is.null(max_ratio)) {
        1.2 * 0.0107 * mz[i] * abs(charge) / 12
      } else max_ratio
      r <- .profile_cor(prof$profiles, i, j)
      ok <- is.finite(ratio) && ratio > 0 && ratio <= bound &&
        (is.null(prof$profiles) || (is.finite(r) && r >= min_cor))
      err <- (mz[j] - mz[i] - spacing) / mz[j] * 1e6
      if (ok && (is.null(best) || abs(err) < abs(best$err))) {
        best <- list(i = i, err = err, ratio = ratio, bound = bound, r = r)
      }
    }
    if (!is.null(best)) {
      out$isotope[j] <- TRUE
      out$parent[j] <- best$i
      out$parent_mz[j] <- mz[best$i]
      out$ppm_error[j] <- best$err
      out$intensity_ratio[j] <- best$ratio
      out$max_ratio[j] <- best$bound
      out$correlation[j] <- best$r
    }
  }
  out
}

.adduct_relations <- function(polarity = c("positive", "negative")) {
  polarity <- match.arg(polarity)
  m <- .spamtp_mass_constants()
  if (polarity == "positive") {
    data.frame(
      relation = c("Na-H", "K-H", "K-Na", "H2O"),
      delta = c(m$sodium_ion - m$proton, m$potassium_ion - m$proton,
                m$potassium_ion - m$sodium_ion, m$water),
      stringsAsFactors = FALSE)
  } else {
    data.frame(relation = c("Cl-H", "H2O"),
               delta = c(m$chloride + m$proton, m$water),
               stringsAsFactors = FALSE)
  }
}

#' Group co-localised peaks that differ by adduct spacings
#'
#' Links two peaks when their m/z difference matches an adduct exchange
#' (for example `[M+Na]+` minus `[M+H]+`) or a water loss within `ppm`, and
#' their spatial profiles correlate. Split peaks of one ion (features closer
#' than `split_ppm` with correlated profiles) are linked as well. Connected
#' components are candidate molecules observed through several ions. Matrix
#' cluster and salt families also form such chains, so peaks flagged by
#' [flagMatrixPeaks()] or [flagIsotopePeaks()] should be passed as `exclude`.
#'
#' @param mz Numeric vector of observed m/z values.
#' @param intensity Features x observations matrix (rows in the order of
#'   `mz`), for example per-pixel normalised intensities.
#' @param ppm Spacing tolerance in parts per million of the heavier peak.
#' @param min_cor Smallest profile correlation for a link.
#' @param polarity Ion polarity; selects the spacing set.
#' @param split_ppm Largest separation treated as a split peak (0 disables).
#' @param exclude Optional logical vector (one value per peak) or peak
#'   indices to leave unlinked.
#'
#' @return A list with `edges` (linked peak pairs, relation, spacing error,
#'   correlation) and `groups` (one row per peak with its group id and size;
#'   singletons, including excluded peaks, get their own group).
#' @examples
#' set.seed(1)
#' base <- abs(rnorm(50))
#' x <- rbind(base, 0.6 * base + rnorm(50, sd = 0.05), abs(rnorm(50)))
#' groupAdductPeaks(c(149.0444, 171.0264, 300.1), x)$groups
#' @export
groupAdductPeaks <- function(mz, intensity, ppm = 3, min_cor = 0.5,
                             polarity = c("positive", "negative"),
                             split_ppm = 0, exclude = NULL) {
  mz <- suppressWarnings(as.numeric(mz))
  if (any(!is.finite(mz) | mz <= 0)) {
    stop("mz must contain positive finite values.", call. = FALSE)
  }
  if (is.null(dim(intensity)) || nrow(intensity) != length(mz)) {
    stop("intensity must be a features x observations matrix.",
         call. = FALSE)
  }
  profiles <- .feature_profiles(mz, intensity)$profiles
  skip <- .exclusion_mask(exclude, length(mz))
  rel <- .adduct_relations(match.arg(polarity))
  ord <- order(mz)
  smz <- mz[ord]
  edges <- list()
  link <- function(i, j, relation, err) {
    if (skip[i] || skip[j]) return(invisible(NULL))
    r <- .profile_cor(profiles, i, j)
    if (is.finite(r) && r >= min_cor) {
      edges[[length(edges) + 1L]] <<- data.frame(
        from = i, to = j, from_mz = mz[i], to_mz = mz[j],
        relation = relation, ppm_error = err, correlation = r,
        stringsAsFactors = FALSE)
    }
  }
  for (i in seq_along(mz)) {
    for (k in seq_len(nrow(rel))) {
      target <- mz[i] + rel$delta[k]
      tol <- target * ppm * 1e-6
      lo <- findInterval(target - tol, smz) + 1L
      hi <- findInterval(target + tol, smz)
      if (hi >= lo) for (j in ord[lo:hi]) {
        link(i, j, rel$relation[k], (mz[j] - target) / target * 1e6)
      }
    }
    if (split_ppm > 0) {
      tol <- mz[i] * split_ppm * 1e-6
      hi <- findInterval(mz[i] + tol, smz)
      lo <- findInterval(mz[i], smz) + 1L
      if (hi >= lo) for (j in ord[lo:hi]) {
        if (j != i && mz[j] > mz[i]) {
          link(i, j, "split", (mz[j] - mz[i]) / mz[i] * 1e6)
        }
      }
    }
  }
  edges <- if (length(edges)) do.call(rbind, edges) else data.frame(
    from = integer(), to = integer(), from_mz = numeric(), to_mz = numeric(),
    relation = character(), ppm_error = numeric(), correlation = numeric())

  # connected components by union-find
  parent <- seq_along(mz)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  for (e in seq_len(nrow(edges))) {
    a <- find(edges$from[e])
    b <- find(edges$to[e])
    if (a != b) parent[max(a, b)] <- min(a, b)
  }
  root <- vapply(seq_along(mz), find, integer(1))
  group <- match(root, unique(root))
  list(edges = edges,
       groups = data.frame(peak = seq_along(mz), mz = mz, group = group,
                           group_size = as.integer(table(group)[group])))
}

.match_mz <- function(x, reference, ppm = 0.5) {
  ord <- order(reference)
  sref <- reference[ord]
  pos <- findInterval(x, sref)
  lower <- pmax(pos, 1L)
  upper <- pmin(pos + 1L, length(sref))
  pick <- ifelse(abs(sref[lower] - x) <= abs(sref[upper] - x), lower, upper)
  ok <- abs(sref[pick] - x) / x * 1e6 <= ppm
  ifelse(ok, ord[pick], NA_integer_)
}

.empty_group_tables <- function() {
  list(
    compounds = data.frame(group = integer(), ramp_id = character(),
                           n_peaks = integer(), n_adducts = integer(),
                           peaks = character(), adducts = character(),
                           score_sum = numeric(), score_max = numeric(),
                           stringsAsFactors = FALSE),
    peaks = data.frame(peak = integer(), mz = numeric(), group = integer(),
                       ramp_ids = character(), n_identities = integer(),
                       support_peaks = integer(), score_sum = numeric(),
                       stringsAsFactors = FALSE))
}

#' Assign identities to co-localised adduct groups
#'
#' Within each group from [groupAdductPeaks()], keeps compounds that explain
#' at least `min_peaks` member peaks through at least two different adducts
#' in a per-peak candidate table, so compounds supported by a single ion are
#' dropped. Each peak then takes the supported compounds that explain the
#' most peaks of its group; several molecules chained into one group keep
#' their own identities, and remaining ties are reported as ambiguity.
#'
#' @param candidates Per-peak candidates with `observed_mz`, `ramp_ids`,
#'   `adduct` and `score` columns, for example from
#'   [queryMZAnnotationIndex()] with `min_score = 0`.
#' @param groups The `groups` element of [groupAdductPeaks()].
#' @param min_peaks Smallest number of member peaks a compound must explain.
#' @param min_score Smallest per-peak score counted as support.
#' @param match_ppm Tolerance for matching candidate m/z to group peaks.
#'
#' @return A list with `compounds` (one row per group and supported compound:
#'   supporting peaks and adducts, summed and maximum peak scores) and
#'   `peaks` (one row per supported peak: its identities, their number and
#'   the number of group peaks they explain).
#' @examples
#' groups <- data.frame(peak = 1:2, mz = c(149.0444, 171.0264), group = 1,
#'                      group_size = 2)
#' cand <- data.frame(observed_mz = c(149.0444, 171.0264),
#'                    ramp_ids = "RAMP_C_1", adduct = c("M+H", "M+Na"),
#'                    score = c(0.4, 0.3))
#' annotateAdductGroups(cand, groups)$peaks
#' @export
annotateAdductGroups <- function(candidates, groups, min_peaks = 2L,
                                 min_score = 0, match_ppm = 0.5) {
  need <- c("observed_mz", "ramp_ids", "adduct", "score")
  if (!all(need %in% names(candidates))) {
    stop("candidates needs columns: ", paste(need, collapse = ", "),
         call. = FALSE)
  }
  if (!all(c("peak", "mz", "group", "group_size") %in% names(groups))) {
    stop("groups must be the groups element of groupAdductPeaks().",
         call. = FALSE)
  }
  min_peaks <- max(2L, as.integer(min_peaks))
  multi <- groups[groups$group_size >= min_peaks, , drop = FALSE]
  cand <- candidates[candidates$score >= min_score &
                       !is.na(candidates$ramp_ids), , drop = FALSE]
  if (!nrow(multi) || !nrow(cand)) return(.empty_group_tables())
  cand$peak <- multi$peak[.match_mz(cand$observed_mz, multi$mz, match_ppm)]
  cand <- cand[!is.na(cand$peak), , drop = FALSE]
  if (!nrow(cand)) return(.empty_group_tables())
  ids <- strsplit(as.character(cand$ramp_ids), ";[[:space:]]*")
  long <- data.frame(peak = rep(cand$peak, lengths(ids)),
                     ramp_id = unlist(ids, use.names = FALSE),
                     adduct = rep(cand$adduct, lengths(ids)),
                     score = rep(cand$score, lengths(ids)),
                     stringsAsFactors = FALSE)
  long <- long[!is.na(long$ramp_id) & nzchar(long$ramp_id), , drop = FALSE]
  long$group <- multi$group[match(long$peak, multi$peak)]
  # best support per (group, compound, peak)
  long <- long[order(-long$score), , drop = FALSE]
  long <- long[!duplicated(long[, c("group", "ramp_id", "peak")]), ,
               drop = FALSE]
  if (!nrow(long)) return(.empty_group_tables())
  key <- paste(long$group, long$ramp_id, sep = "\r")
  compounds <- do.call(rbind, lapply(split(seq_len(nrow(long)), key),
                                     function(i) {
    i <- i[order(long$peak[i])]
    data.frame(group = long$group[i[1]], ramp_id = long$ramp_id[i[1]],
               n_peaks = length(unique(long$peak[i])),
               n_adducts = length(unique(long$adduct[i])),
               peaks = paste(long$peak[i], collapse = ";"),
               adducts = paste(long$adduct[i], collapse = ";"),
               score_sum = sum(long$score[i]), score_max = max(long$score[i]),
               stringsAsFactors = FALSE)
  }))
  compounds <- compounds[compounds$n_peaks >= min_peaks &
                           compounds$n_adducts >= 2, , drop = FALSE]
  if (!nrow(compounds)) return(.empty_group_tables())
  compounds <- compounds[order(compounds$group, -compounds$n_peaks,
                               -compounds$score_sum), , drop = FALSE]
  rownames(compounds) <- NULL

  # per peak: the supported compounds explaining the most group peaks
  members <- strsplit(compounds$peaks, ";")
  per_peak <- data.frame(
    peak = as.integer(unlist(members, use.names = FALSE)),
    row = rep(seq_len(nrow(compounds)), lengths(members)))
  peaks <- do.call(rbind, lapply(split(per_peak$row, per_peak$peak),
                                 function(rows) {
    support <- compounds$n_peaks[rows]
    best <- rows[support == max(support)]
    data.frame(ramp_ids = paste(sort(unique(compounds$ramp_id[best])),
                                collapse = "; "),
               n_identities = length(unique(compounds$ramp_id[best])),
               support_peaks = max(support),
               score_sum = max(compounds$score_sum[best]),
               stringsAsFactors = FALSE)
  }))
  peaks$peak <- as.integer(rownames(peaks))
  peaks$mz <- groups$mz[match(peaks$peak, groups$peak)]
  peaks$group <- groups$group[match(peaks$peak, groups$peak)]
  peaks <- peaks[order(peaks$group, peaks$peak),
                 c("peak", "mz", "group", "ramp_ids", "n_identities",
                   "support_peaks", "score_sum")]
  rownames(peaks) <- NULL
  list(compounds = compounds, peaks = peaks)
}

#' Restrict a chemical table to pathway-connected compounds
#'
#' Subsets a RaMP chemical table to compounds linked to pathways from the
#' chosen sources, for use with [buildMZAnnotationIndex()]. A smaller search
#' space can reduce chance matches while losing relevant compounds, including
#' lipid species absent from those pathways. Compare coverage and mass-shift
#' diagnostics with [estimateAnnotationFDR()], and justify the database scope
#' independently of the results. Pathway membership is not evidence that a
#' compound occurs in the sample; selecting references by pathway membership
#' also changes the background of downstream pathway analyses.
#'
#' @param chem_props RaMP chemical table with a `ramp_id` column.
#' @param analytehaspathway RaMP analyte-pathway table with `rampId` and
#'   `pathwaySource` columns.
#' @param sources Pathway sources to accept; the default keeps manually
#'   curated databases.
#'
#' @return The subset of `chem_props`.
#' @examples
#' chem <- data.frame(ramp_id = c("RAMP_C_1", "RAMP_C_2"))
#' ahp <- data.frame(rampId = "RAMP_C_1", pathwaySource = "reactome")
#' annotationDatabaseScope(chem, ahp)
#' @export
annotationDatabaseScope <- function(chem_props, analytehaspathway,
                                    sources = c("reactome", "kegg", "wiki")) {
  if (!"ramp_id" %in% names(chem_props)) {
    stop("chem_props needs a ramp_id column.", call. = FALSE)
  }
  if (!all(c("rampId", "pathwaySource") %in% names(analytehaspathway))) {
    stop("analytehaspathway needs rampId and pathwaySource columns.",
         call. = FALSE)
  }
  keep <- unique(analytehaspathway$rampId[
    analytehaspathway$pathwaySource %in% sources])
  chem_props[chem_props$ramp_id %in% keep, , drop = FALSE]
}
