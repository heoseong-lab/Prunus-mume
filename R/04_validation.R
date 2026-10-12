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
temporal <- bind_rows(lapply(sort(unique(dat$year)), function(yy) {
  a <- profile_and_weights(dat[dat$year != yy, ]); test <- dat[dat$year == yy, ]
  common <- sort(intersect(rownames(a$x), test$no)); test <- test[match(common, test$no), ]
  xx <- data.frame(row.names = common)
  for (tr in CONT) xx[[tr]] <- test[[paste0(tr, "_adj")]]
  for (tr in NOM)  xx[[tr]] <- factor(test[[tr]])
  for (tr in ORD)  xx[[tr]] <- ord_score(test[[tr]], tr)
  ws <- scheme_weights(a)
  bind_rows(lapply(names(ws), function(s) {
    tr_run <- cluster_run(a$x[common, ], ws[[s]], k = K); te_run <- cluster_run(xx, ws[[s]], k = K)
    data.frame(held_out_year = yy, scheme = s, n_common = length(common), ARI = ari(tr_run$cl, te_run$cl),
               distance_spearman = cor(as.vector(tr_run$gower), as.vector(te_run$gower), method = "spearman"))
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
