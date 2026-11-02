# =============================================================================
# 04_validation.R  Disjoint splits (weight transfer) and leave-one-year-out
# Three weighting schemes are carried through every check:
#   weighted  = primary weights (between-year correlation / kappa)
#   equal     = equal weights for the measured variables (reference-assessed descriptors excluded)
#   retained  = equal weights for the variables with positive primary weight
# =============================================================================
MEAS <- setdiff(PRIMARY, FIXED_TRAITS)   # variables with a measured between-year coefficient
scheme_weights <- function(a) list(weighted = a$weights[PRIMARY], equal = setNames(as.numeric(!(PRIMARY %in% FIXED_TRAITS)), PRIMARY),
                                   retained = setNames(as.numeric(a$weights[PRIMARY] > 0), PRIMARY))
ref_part <- list(weighted = base_run$cl, equal = cutree(runs$Equal_all$hc, K), retained = cutree(runs$Equal_retained$hc, K))

errors <- list()

# ---- disjoint splits -----------------------------------------------------------
ids <- sort(unique(dat$no)); set.seed(SEED + 2); split_rows <- list()
for (b in seq_len(B_SPLIT)) {
  train <- sample(ids, floor(length(ids) / 2)); test <- setdiff(ids, train)
  attempt <- try({
    a <- profile_and_weights(dat[dat$no %in% train, ]); v <- profile_and_weights(dat[dat$no %in% test, ])
    trans <- cluster_run(v$x, a$weights[PRIMARY], k = K); local <- cluster_run(v$x, v$weights[PRIMARY], k = K)
    eq <- cluster_run(v$x, setNames(as.numeric(!(PRIMARY %in% FIXED_TRAITS)), PRIMARY), k = K)
    rt <- cluster_run(v$x, setNames(as.numeric(a$weights[PRIMARY] > 0), PRIMARY), k = K)
    split_rows[[b]] <- data.frame(iteration = b, n_train = length(train), n_test = length(test),
      ARI_transfer_local = ari(trans$cl, local$cl), ARI_transfer_equal = ari(trans$cl, eq$cl), ARI_transfer_retained = ari(trans$cl, rt$cl),
      weight_spearman = cor(a$weights[MEAS], v$weights[MEAS], method = "spearman"))
    TRUE
  }, silent = TRUE)
  if (inherits(attempt, "try-error")) errors[[length(errors) + 1]] <- data.frame(analysis = "split", iteration = b, error = as.character(attempt))
}
split_result <- bind_rows(split_rows); save_table(split_result, "disjoint_weight_transfer")

# ---- leave-one-year-out --------------------------------------------------------
# Training: profiles and weights from two years. Test: the held-out year's stage-1 adjusted records.
precompute_components <- function(x) {
  n <- nrow(x); dj_list <- list(); ok_list <- list()
  for (tr in PRIMARY) {
    v <- x[[tr]]; ok <- outer(!is.na(v), !is.na(v), "&")
    if (is.factor(v)) dj <- outer(as.character(v), as.character(v), "!=") * 1
    else { ran <- diff(range(v, na.rm = TRUE)); dj <- if (ran > 0) abs(outer(v, v, "-")) / ran else matrix(0, n, n) }
    dj[!ok] <- 0; dj_list[[tr]] <- dj[upper.tri(dj)]; ok_list[[tr]] <- (ok * 1)[upper.tri(ok)]
  }
  list(dj = do.call(cbind, dj_list), ok = do.call(cbind, ok_list), n = n)
}
fast_gower_vec <- function(comp, w) {
  w <- w[colnames(comp$dj)]; keep <- which(is.finite(w) & w > 0)
  as.numeric((comp$dj[, keep, drop = FALSE] %*% w[keep]) / pmax(comp$ok[, keep, drop = FALSE] %*% w[keep], 1e-15))
}

set.seed(SEED + 10)
temporal <- bind_rows(lapply(sort(unique(dat$year)), function(yy) {
  a <- profile_and_weights(dat[dat$year != yy, ]); test <- dat[dat$year == yy, ]
  common <- sort(intersect(rownames(a$x), test$no)); test <- test[match(common, test$no), ]
  xx <- data.frame(row.names = common)
  for (tr in CONT) xx[[tr]] <- test[[paste0(tr, "_adj")]]
  for (tr in NOM)  xx[[tr]] <- factor(test[[tr]])
  for (tr in ORD)  xx[[tr]] <- ord_score(test[[tr]], tr)
  ws <- scheme_weights(a)
  runs_s <- lapply(names(ws), function(s) {
    list(tr = cluster_run(a$x[common, ], ws[[s]], k = K), te = cluster_run(xx, ws[[s]], k = K))
  })
  names(runs_s) <- names(ws)
  m <- length(common); ut_base <- which(upper.tri(matrix(0, m, m)), arr.ind = TRUE); r1 <- ut_base[, 1]; r2 <- ut_base[, 2]
  D_mats <- lapply(runs_s, function(z) list(tr = as.matrix(z$tr$gower), te = as.matrix(z$te$gower)))
  boot_mat <- replicate(B_BOOT, {
    idx <- sample.int(m, m, replace = TRUE); i1 <- idx[r1]; i2 <- idx[r2]; kp <- i1 != i2
    lin <- (i2[kp] - 1L) * m + i1[kp]
    vapply(names(ws), function(s) cor(D_mats[[s]]$tr[lin], D_mats[[s]]$te[lin], method = "spearman"), numeric(1))
  })
  comp_tr <- precompute_components(a$x[common, ]); comp_te <- precompute_components(xx)
  pos_tr  <- names(ws$weighted)[ws$weighted > 0]
  perm_rhos <- replicate(B_BOOT, {
    wp <- ws$weighted; wp[pos_tr] <- sample(ws$weighted[pos_tr])
    cor(fast_gower_vec(comp_tr, wp), fast_gower_vec(comp_te, wp), method = "spearman")
  })
  obs_w_rho <- cor(as.vector(runs_s$weighted$tr$gower), as.vector(runs_s$weighted$te$gower), method = "spearman")
  bind_rows(lapply(names(ws), function(s) {
    tr_run <- runs_s[[s]]$tr; te_run <- runs_s[[s]]$te
    r_obs <- cor(as.vector(tr_run$gower), as.vector(te_run$gower), method = "spearman")
    diff_w <- boot_mat["weighted", ] - boot_mat[s, ]
    data.frame(held_out_year = yy, scheme = s, n_common = m, ARI = ari(tr_run$cl, te_run$cl),
               distance_spearman = r_obs,
               rho_low = quantile(boot_mat[s, ], 0.025, na.rm = TRUE),
               rho_high = quantile(boot_mat[s, ], 0.975, na.rm = TRUE),
               diff_vs_scheme = if (s == "weighted") NA_real_ else obs_w_rho - r_obs,
               diff_low = if (s == "weighted") NA_real_ else quantile(diff_w, 0.025, na.rm = TRUE),
               diff_high = if (s == "weighted") NA_real_ else quantile(diff_w, 0.975, na.rm = TRUE),
               perm_median = if (s == "weighted") median(perm_rhos) else NA_real_,
               perm_low = if (s == "weighted") quantile(perm_rhos, 0.025) else NA_real_,
               perm_high = if (s == "weighted") quantile(perm_rhos, 0.975) else NA_real_,
               perm_p = if (s == "weighted") (sum(perm_rhos >= obs_w_rho) + 1) / (length(perm_rhos) + 1) else NA_real_)
  }))
}))
save_table(temporal, "leave_one_year_out")
save_table(if (length(errors)) bind_rows(errors) else data.frame(analysis = character(), iteration = integer(), error = character()), "validation_failures")
validation_summary <- list(
  splits = list(requested = B_SPLIT, completed = nrow(split_result), transfer_local = q95(split_result$ARI_transfer_local),
    transfer_equal = q95(split_result$ARI_transfer_equal), transfer_retained = q95(split_result$ARI_transfer_retained),
    weight_spearman = q95(split_result$weight_spearman), share_local_gt_equal = mean(split_result$ARI_transfer_local > split_result$ARI_transfer_equal)),
  temporal = temporal)
write_json(validation_summary, file.path(OUT, "validation_summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 4)
stopifnot(nrow(split_result) >= .9 * B_SPLIT)
