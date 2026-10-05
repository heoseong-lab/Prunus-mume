# =============================================================================
# 00_config.R  Packages, constants, trait definitions, shared functions
#
# Two-stage reproducibility-weighted analysis of the P. mume collection.
#   Stage 1 (01_prepare.R): within each year, fruit traits are adjusted by the
#           harvest-batch median and the other quantitative traits by the year
#           median, so that values from different years are on one scale.
#   Stage 2 (02_reproducibility.R): the accession profile is the mean of the
#           adjusted values over the observed years, and the reproducibility of
#           a quantitative trait is the mean between-year Spearman correlation
#           of the adjusted values. Categorical traits use Fleiss' kappa.
#           Both weights are single-observation coefficients, so they are
#           on the same basis.
# Exclusions requested for this analysis: the two disease descriptors and the
# accessions KJM261 and KJM264 (petal number recorded as 3).
# =============================================================================
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(cluster); library(ggplot2); library(dplyr); library(tidyr)
  library(readr); library(jsonlite); library(ordinal)
})

SEED     <- 20261001L
B_BOOT   <- as.integer(Sys.getenv("PM_BOOT",   "1000"))  # accession bootstraps for the reproducibility coefficients
B_SPLIT  <- as.integer(Sys.getenv("PM_SPLIT",  "100"))   # disjoint splits
B_RANDOM <- as.integer(Sys.getenv("PM_RANDOM", "1000"))  # random benchmark sets
set.seed(SEED)

EXCLUDE_ACCESSIONS <- c("KJM261", "KJM264")
EXCLUDE_TRAITS     <- c("brown_rot_severity", "shot_hole_severity")
# Reference-year harmonization: descriptors taken from a single reference assessment per accession.
# For each descriptor, the first listed year in which the accession has a record supplies the value
# for all of its other years (fruit shape: 2024, otherwise the earliest year; flesh color: 2026,
# otherwise the latest year). The harmonized file is written to data/phenotypes_harmonized.csv
# together with the reference year used for each record (columns <trait>_reference_year).
HARMONIZE    <- list(fruit_shape = c(2024L, 2025L, 2026L), flesh_color = c(2026L, 2025L, 2024L))
HARMONIZE_FALLBACK <- Sys.getenv("PM_FALLBACK", "nearest")   # "nearest": use the nearest listed year when the first is absent; "missing": set the descriptor to NA
# The harmonized descriptors are identical across years by construction, so they have no measured
# between-year coefficient. They are used descriptively (cluster profiles, Table S2, priority-set
# coverage) and enter the distance with the weight given here (0 = not part of the distance).
FIXED_WEIGHTS <- c(fruit_shape = as.numeric(Sys.getenv("PM_W_FRUIT_SHAPE", "0")), flesh_color = as.numeric(Sys.getenv("PM_W_FLESH_COLOR", "0")))
FIXED_TRAITS  <- names(FIXED_WEIGHTS)

OUT <- file.path(ROOT, "output")
for (p in c("tables", "figures", "objects")) dir.create(file.path(OUT, p), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) { write_csv(as.data.frame(x), file.path(OUT, "tables", paste0(name, ".csv")), na = ""); invisible(x) }
save_plot  <- function(p, name, w = 7.3, h = 5) {
  ggsave(file.path(OUT, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = 400, bg = "white")
  ggsave(file.path(OUT, "figures", paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
}
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom"))
# ---- ShinyCore (Kim et al. 2023, Plant Methods 19:106; https://github.com/heoseong/ShinyCore) ------------
# Port of the selection algorithm of ShinyCore.R (covering phase + thickening phase) to a plain
# function. x: data.frame of samples (rows) x markers (columns) with categorical classes; NA allowed.
#   covering  : start with the sample with the fewest missing values; then repeatedly add the sample
#               that carries the largest number of not-yet-covered classes (classes already covered are
#               masked as NA, markers fully covered are dropped) until the coverage CV reaches cov_max.
#   thickening: keep the first N0 samples of the covering sequence (CV > cov_min) and fill the
#               remaining N1 - N0 places with the samples of highest rarity score
#               RS = sum over markers of (1 - p)/p, p = frequency of the sample's class.
# Returns the core, the covering history and the ShinyCore report values (CV, SH, RS).
CORE_COV_MAX <- as.numeric(Sys.getenv("PM_COV_MAX", "0.99"))   # ShinyCore default: maximum coverage 99%
CORE_COV_MIN <- as.numeric(Sys.getenv("PM_COV_MIN", "0.98"))   # ShinyCore default: minimum coverage 98%
N_CLASSES    <- as.integer(Sys.getenv("PM_CLASSES", "5"))      # equal-width classes for each quantitative profile
sc_count   <- function(v) length(unique(v[!is.na(v)]))
sc_shannon <- function(v) { p <- as.numeric(prop.table(table(v))); -sum(log(p^p, max(length(p), 2))) }
sc_odds    <- function(v) { p <- prop.table(table(v)); (1 - p) / p }
shinycore <- function(x, cov_max = CORE_COV_MAX, cov_min = CORE_COV_MIN) {
  x <- as.data.frame(lapply(x, as.character), row.names = rownames(x), stringsAsFactors = FALSE)
  x <- x[, vapply(x, sc_count, integer(1)) >= 2, drop = FALSE]                       # markers with >= 2 classes
  x <- x[, order(vapply(x, function(v) min(prop.table(table(v))), numeric(1))), drop = FALSE]   # rare markers first
  m <- as.matrix(x); k <- ncol(m); count_all <- apply(m, 2, sc_count)
  noncore <- m; core <- m[0, , drop = FALSE]; history <- character(0); cov_hist <- numeric(0)
  cur_all <- count_all
  repeat {
    i <- which.min(rowSums(is.na(noncore)))                   # fewest NA = most not-yet-covered classes
    new <- noncore[i, , drop = FALSE]; core <- rbind(core, new); noncore <- noncore[-i, , drop = FALSE]
    history <- c(history, rownames(new))
    count_core <- apply(core, 2, sc_count)
    cov_hist <- c(cov_hist, (sum(count_core / cur_all) + (k - length(cur_all))) / k)
    done <- which(count_core == cur_all)
    if (length(done)) { core <- core[, -done, drop = FALSE]; noncore <- noncore[, -done, drop = FALSE]; new <- new[, -done, drop = FALSE]; cur_all <- cur_all[-done] }
    if (tail(cov_hist, 1) >= cov_max || nrow(noncore) == 0 || ncol(noncore) == 0) break
    hit <- !is.na(noncore) & noncore == matrix(new, nrow(noncore), ncol(noncore), byrow = TRUE)   # mask covered classes
    noncore[hit] <- NA
  }
  N1 <- length(history); N0 <- min(which(cov_hist > cov_min))
  odds <- lapply(as.data.frame(m, stringsAsFactors = FALSE), sc_odds)
  score <- rowSums(sapply(seq_len(k), function(j) { o <- odds[[j]]; v <- as.numeric(o[match(m[, j], names(o))]); ifelse(is.na(v), 0, v) }))
  names(score) <- rownames(m); score <- sort(score, decreasing = TRUE)
  core_ids <- c(history[seq_len(N0)], names(score)[!names(score) %in% history[seq_len(N0)]][seq_len(N1 - N0)])
  sub <- m[core_ids, , drop = FALSE]
  list(core = core_ids, history = history, coverage_history = cov_hist, N0 = N0, N1 = N1, score = score,
       CV = mean(apply(sub, 2, sc_count) / count_all), SH = mean(apply(sub, 2, sc_shannon)), RS = log10(sum(score[core_ids])),
       SH_all = mean(apply(m, 2, sc_shannon)), n_markers = k, n_classes = sum(count_all))
}
# equal-width classes of a quantitative profile over its observed range (class labels 1..n)
sc_classes <- function(v, n = N_CLASSES) {
  r <- range(v, na.rm = TRUE); if (diff(r) == 0) return(rep("1", length(v)))
  cl <- pmin(floor((v - r[1]) / diff(r) * n) + 1, n); ifelse(is.na(cl), NA_character_, as.character(cl))
}

CL_COL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#332288", "#888888", "#882255", "#44AA99")

# ---- trait sets -------------------------------------------------------------
CONT       <- c("fruit_weight", "petal_number", "soluble_solids_content", "titratable_acidity", "bloom_doy")
FRUIT_CONT <- c("fruit_weight", "soluble_solids_content", "titratable_acidity")   # adjusted by harvest batch
NOM        <- c("fruit_shape", "ground_color", "flesh_color", "leaf_base_shape")
LEV <- list(
  skin_color_intensity = c("absent or very weak", "weak", "medium", "strong", "very strong"),
  stone_adherence      = c("clingstone", "semi-clingstone", "semi-freestone", "freestone"),
  leaf_green_intensity = c("light", "medium", "dark"),
  leaf_pubescence      = c("absent", "weak", "medium", "strong"),
  leaf_angle           = c("acute", "right-angled", "moderately obtuse", "strongly obtuse"),
  leaf_length          = c("short", "medium", "long"),
  leaf_width           = c("narrow", "medium", "broad"),
  flower_diameter      = c("small", "medium", "large"))
ORD     <- names(LEV)
CAT     <- c(NOM, ORD)
PRIMARY <- c(CONT, CAT)                       # 17 variables
LABEL <- c(fruit_weight = "Fruit weight", petal_number = "Petal number", soluble_solids_content = "Soluble solids content",
           titratable_acidity = "Titratable acidity", bloom_doy = "Full bloom date", fruit_shape = "Fruit shape",
           ground_color = "Ground color", flesh_color = "Flesh color", leaf_base_shape = "Leaf base shape",
           skin_color_intensity = "Skin color intensity", stone_adherence = "Stone adherence",
           leaf_green_intensity = "Leaf green intensity", leaf_pubescence = "Leaf pubescence", leaf_angle = "Leaf angle",
           leaf_length = "Leaf length", leaf_width = "Leaf width", flower_diameter = "Flower diameter")
UNIT  <- c(fruit_weight = "g", petal_number = "petals per flower", soluble_solids_content = "°Brix",
           titratable_acidity = "%", bloom_doy = "day of year")
DOUBLE_MIN <- 10   # mean petal number at or above which an accession is called double-flowered

# ---- small helpers ----------------------------------------------------------
mode_value <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_character_); tb <- table(x); sort(names(tb)[tb == max(tb)])[1] }
mode_tie   <- function(x) { tb <- table(x); length(tb) > 1 && sum(tb == max(tb)) > 1 }
ord_score  <- function(v, tr) 1 + 8 * (match(v, LEV[[tr]]) - 1) / (length(LEV[[tr]]) - 1)   # equally spaced 1-9
q95 <- function(x) quantile(x, c(.025, .5, .975), na.rm = TRUE, names = FALSE)

# Adjusted Rand index (Hubert & Arabie 1985)
ari <- function(a, b) {
  tab <- table(a, b); ch <- function(x) x * (x - 1) / 2; n <- sum(tab)
  if (n < 2) return(NA_real_)
  A <- sum(ch(rowSums(tab))); B <- sum(ch(colSums(tab))); E <- A * B / ch(n); den <- (A + B) / 2 - E
  if (abs(den) < 1e-12) return(ifelse(identical(outer(a, a, "=="), outer(b, b, "==")), 1, 0))
  (sum(ch(tab)) - E) / den
}

# Fleiss' kappa with years as raters; uses accessions with the maximum common number of ratings
fleiss <- function(id, x) {
  ok <- !is.na(x); id <- as.character(id[ok]); x <- as.character(x[ok]); n <- table(id)
  if (!length(n)) return(c(value = NA, n = 0, raters = 0))
  m <- max(n); ids <- names(n)[n == m]
  if (m < 2) return(c(value = NA, n = length(ids), raters = m))
  tab <- table(factor(id[id %in% ids], levels = ids), x[id %in% ids]); pe <- sum((colSums(tab) / sum(tab))^2)
  k <- if (pe >= 1 - 1e-12) NA_real_ else (mean((rowSums(tab^2) - m) / (m * (m - 1))) - pe) / (1 - pe)
  c(value = k, n = nrow(tab), raters = m)
}

# Krippendorff's alpha (coincidence matrix); ordinal metric for ordinal descriptors
alpha_df <- function(d, tr) {
  a <- data.frame(id = as.character(d$no), x = as.character(d[[tr]])); a <- a[!is.na(a$x), ]
  ns <- table(a$id); a <- a[a$id %in% names(ns)[ns >= 2], ]
  lv <- if (tr %in% ORD) LEV[[tr]] else sort(unique(a$x)); q <- length(lv)
  if (q < 2 || !nrow(a)) return(NA_real_)
  O <- matrix(0, q, q)
  for (v in split(a$x, a$id)) { n <- tabulate(match(v, lv), q); m <- sum(n); if (m > 1) O <- O + (outer(n, n) - diag(n)) / (m - 1) }
  marg <- rowSums(O); nt <- sum(marg); if (nt <= 1) return(NA_real_)
  delta <- 1 - diag(q)
  if (tr %in% ORD) for (i in seq_len(q)) for (j in seq_len(q)) if (i != j) { lo <- min(i, j); hi <- max(i, j); delta[i, j] <- (sum(marg[lo:hi]) - (marg[lo] + marg[hi]) / 2)^2 }
  E <- (outer(marg, marg) - diag(marg)) / (nt - 1); den <- sum(E * delta)
  if (den <= 0) NA_real_ else 1 - sum(O * delta) / den
}

# Between-year Spearman correlation of a (stage-1 adjusted) quantitative trait.
# Returns one row per year pair and the mean over pairs with at least MIN_PAIR common accessions.
MIN_PAIR <- 20
year_pair_cor <- function(d, tr, col = tr) {
  yrs <- sort(unique(d$year)); if (length(yrs) < 2) return(NULL)
  prs <- combn(yrs, 2, simplify = FALSE)
  bind_rows(lapply(prs, function(yy) {
    a <- d[d$year == yy[1], c("no", col)]; b <- d[d$year == yy[2], c("no", col)]
    z <- merge(a, b, by = "no"); z <- z[complete.cases(z), ]
    data.frame(trait = tr, year1 = yy[1], year2 = yy[2], n = nrow(z),
               rho = if (nrow(z) >= MIN_PAIR) cor(z[, 2], z[, 3], method = "spearman") else NA_real_)
  }))
}
mean_rho <- function(pc) if (is.null(pc) || all(is.na(pc$rho))) NA_real_ else mean(pc$rho, na.rm = TRUE)

# Latent-scale intraclass correlation from a cumulative link mixed model with year as a fixed effect
latent_icc <- function(d, tr) {
  a <- d[!is.na(d[[tr]]), ]; ids <- names(which(table(a$no) >= 2)); a <- a[a$no %in% ids, ]
  if (length(unique(a$year)) < 2) return(NA_real_)
  a$y <- ordered(a[[tr]], levels = LEV[[tr]]); a$no <- factor(a$no); a$year_f <- factor(a$year)
  fit <- try(suppressWarnings(ordinal::clmm(y ~ year_f + (1 | no), data = a, Hess = FALSE)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  v <- as.numeric(ordinal::VarCorr(fit)$no[1]); v / (v + pi^2 / 3)
}

sb <- function(r, n) n * r / (1 + (n - 1) * r)   # Spearman-Brown step-up (sensitivity analysis only)
