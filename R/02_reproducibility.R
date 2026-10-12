# =============================================================================
# 02_reproducibility.R  Stage 2: accession profiles and reproducibility weights
#
# profile_and_weights(d) is the single function reused by every resampling
# analysis. It takes accession-year records (already stage-1 adjusted) and
# returns the accession profile matrix, the weights, and a reliability table.
#   quantitative trait : profile = mean of adjusted values over observed years
#                        reproducibility = mean between-year Spearman correlation
#   nominal trait      : profile = modal level; reproducibility = Fleiss' kappa
#   ordinal trait      : profile = median of 1-9 scores; reproducibility = Fleiss' kappa
# Weight = max(coefficient, 0). No step-up is applied to either data type.
# =============================================================================
profile_and_weights <- function(d, latent = FALSE) {
  ids <- sort(unique(d$no))
  x <- data.frame(row.names = ids)
  rel <- list()
  for (tr in CONT) {
    col <- paste0(tr, "_adj")
    m <- as.numeric(tapply(d[[col]], factor(d$no, levels = ids), mean, na.rm = TRUE)[ids]); m[is.nan(m)] <- NA
    x[[tr]] <- m
    pc <- year_pair_cor(d, tr, col)
    r  <- mean_rho(pc)
    rel[[tr]] <- data.frame(trait = tr, type = "quantitative", fixed = FALSE, n_records = sum(!is.na(d[[col]])),
                            n_accessions = sum(!is.na(x[[tr]])), coefficient = r, alpha = NA, latent_icc = NA,
                            weight = ifelse(is.finite(r), max(0, r), 0))
  }
  for (tr in CAT) {
    f  <- fleiss(d$no, d[[tr]])
    if (tr %in% NOM) {
      x[[tr]] <- factor(vapply(ids, function(i) mode_value(d[[tr]][d$no == i]), character(1)), levels = sort(unique(na.omit(d[[tr]]))))
    } else {
      sc <- ord_score(d[[tr]], tr)
      x[[tr]] <- vapply(ids, function(i) if (all(is.na(sc[d$no == i]))) NA_real_ else median(sc[d$no == i], na.rm = TRUE), numeric(1))
    }
    k <- f[["value"]]; fx <- tr %in% FIXED_TRAITS
    rel[[tr]] <- data.frame(trait = tr, type = ifelse(tr %in% NOM, "nominal", "ordinal"), fixed = fx, n_records = sum(!is.na(d[[tr]])),
                            n_accessions = sum(!is.na(x[[tr]])), coefficient = if (fx) NA_real_ else k, alpha = if (fx) NA_real_ else alpha_df(d, tr),
                            latent_icc = if (latent && tr %in% ORD) latent_icc(d, tr) else NA,
                            weight = if (fx) unname(FIXED_WEIGHTS[tr]) else ifelse(is.finite(k), max(0, k), 0))
  }
  rel <- bind_rows(rel)
  list(x = x, weights = setNames(rel$weight, rel$trait), reliability = rel)
}

main <- profile_and_weights(dat)
pair_cors <- bind_rows(lapply(CONT, function(tr) year_pair_cor(dat, tr, paste0(tr, "_adj"))))
save_table(pair_cors, "between_year_correlations")
save_table(main$reliability, "trait_reproducibility")
save_table(data.frame(no = rownames(main$x), main$x, check.names = FALSE), "accession_profiles")
save_table(bind_rows(lapply(CAT, function(tr) data.frame(trait = tr, no = unique(dat$no),
  mode_tied = vapply(unique(dat$no), function(i) mode_tie(dat[[tr]][dat$no == i]), logical(1))))), "categorical_mode_ties")

# Bootstrap intervals for the reproducibility coefficients: accessions are drawn with replacement together
# with all their records (duplicates renamed so that each draw is a separate accession), and the between-year
# Spearman correlation (quantitative) or Fleiss' kappa (categorical, non-fixed) is recomputed on each resample.
coefficient_bootstrap <- function(d, B = B_BOOT) {
  ids <- sort(unique(d$no)); by_id <- split(d, d$no); meas <- setdiff(PRIMARY, FIXED_TRAITS)
  bind_rows(lapply(seq_len(B), function(b) {
    draw <- sample(ids, length(ids), replace = TRUE)
    bd <- bind_rows(lapply(seq_along(draw), function(j) { a <- by_id[[draw[j]]]; a$no <- sprintf("B%03d", j); a }))
    data.frame(iteration = b, trait = meas, coefficient = vapply(meas, function(tr)
      if (tr %in% CONT) mean_rho(year_pair_cor(bd, tr, paste0(tr, "_adj"))) else unname(fleiss(bd$no, bd[[tr]])[["value"]]), numeric(1)))
  }))
}
set.seed(SEED + 1)
boot_coef <- coefficient_bootstrap(dat, B_BOOT)
boot_ci <- boot_coef %>% group_by(trait) %>% summarise(n = sum(is.finite(coefficient)),
  coef_low = quantile(coefficient, .025, na.rm = TRUE), coef_high = quantile(coefficient, .975, na.rm = TRUE),
  zero_frequency = mean(pmax(coefficient, 0) == 0, na.rm = TRUE), .groups = "drop")
save_table(boot_coef, "trait_bootstrap_estimates"); save_table(boot_ci, "trait_bootstrap_intervals")

# Year effects on ordinal descriptors (cumulative link mixed models, accessions observed in all three years)
year_shift <- bind_rows(lapply(ORD, function(tr) {
  ids <- names(which(table(dat$no[!is.na(dat[[tr]])]) == 3)); a <- dat[dat$no %in% ids, ]
  a$y <- ordered(a[[tr]], levels = LEV[[tr]]); a$no <- factor(a$no); a$year_f <- factor(a$year)
  full <- try(suppressWarnings(ordinal::clmm(y ~ year_f + (1 | no), data = a, Hess = TRUE)), silent = TRUE)
  null <- try(suppressWarnings(ordinal::clmm(y ~ 1 + (1 | no), data = a, Hess = TRUE)), silent = TRUE)
  p <- if (inherits(full, "try-error") || inherits(null, "try-error")) NA_real_ else
       pchisq(2 * (as.numeric(logLik(full)) - as.numeric(logLik(null))), 2, lower.tail = FALSE)
  codes <- ord_score(a[[tr]], tr); means <- tapply(codes, a$year, mean, na.rm = TRUE)
  data.frame(trait = tr, n_accessions = length(ids), annual_score_range = diff(range(means)), year_LR_p = p)
})) %>% mutate(year_BH_p = p.adjust(year_LR_p, "BH"))
save_table(year_shift, "ordinal_year_effects")

print(main$reliability[, c("trait", "type", "fixed", "n_accessions", "coefficient", "alpha", "weight")] %>%
        left_join(boot_ci[, c("trait", "coef_low", "coef_high")], by = "trait"), digits = 3)
