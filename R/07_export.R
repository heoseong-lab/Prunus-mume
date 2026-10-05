# =============================================================================
# 07_export.R  Numbers quoted in the manuscript, collected in one JSON file
# =============================================================================
vs <- validation_summary; cs <- fromJSON(file.path(OUT, "core_summary.json"))
rel <- main$reliability; ci <- read_csv(file.path(OUT, "tables", "trait_bootstrap_intervals.csv"), show_col_types = FALSE)
nums <- list(
  n_accessions = nrow(main$x), n_records = nrow(dat), records_by_year = as.list(table(dat$year)),
  years_per_accession = as.list(table(years_by_id$n_years)), mean_years = mean(years_by_id$n_years),
  n_boot = B_BOOT, n_primary = length(PRIMARY), n_positive = sum(primary_w > 0), n_zero = sum(primary_w == 0),
  n_double = sum(petal_audit$flower_type == "double"),
  fixed_traits = FIXED_TRAITS, fixed_weights = as.list(FIXED_WEIGHTS), n_harmonized = as.list(table(factor(harm_log$trait, levels = FIXED_TRAITS))),
  reproducibility = setNames(lapply(PRIMARY, function(tr) list(coefficient = rel$coefficient[rel$trait == tr], fixed = rel$fixed[rel$trait == tr], weight = unname(primary_w[tr]),
    alpha = rel$alpha[rel$trait == tr], low = if (tr %in% ci$trait) ci$coef_low[ci$trait == tr] else NA, high = if (tr %in% ci$trait) ci$coef_high[ci$trait == tr] else NA,
    contribution_pct = weight_share$contribution_pct[weight_share$trait == tr], weight_share_pct = weight_share$weight_share_pct[weight_share$trait == tr])), PRIMARY),
  pair_cors = pair_cors,
  ordinal_year_effects = year_shift,
  k = K, silhouette = base_run$sil, silhouette_scan = base_run$scan, cluster_sizes = as.integer(table(base_run$cl)),
  lingoes = base_run$correction, negative_fraction_pct = 100 * base_run$negative_fraction, pcoa_pct = pc_pct[1:2],
  cluster_profiles = profile_c, cluster_categorical = profile_cat, xtab = xtab, xtab_fruit = xtab_fruit, ari_two_way = ari_two_way, ari_cat = ari_cat,
  chi_years = list(statistic = unname(chi_years$statistic), p = chi_years$p.value),
  sensitivity = sens, validation = vs, core = cs, range_tab = range_tab, candidates = candidates,
  no_disease_removed = EXCLUDE_TRAITS, excluded_accessions = EXCLUDE_ACCESSIONS)
write_json(nums, file.path(OUT, "manuscript_numbers.json"), pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 5)
message("manuscript_numbers.json written")
