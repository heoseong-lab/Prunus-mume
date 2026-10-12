# =============================================================================
# 05_priority.R  Core collection selected with the ShinyCore algorithm
#   (Kim, Kim, Moyle & Heo 2023, Plant Methods 19:106; function shinycore() in 00_config.R)
#
# "Markers" for ShinyCore are the accession profiles of the 17 variables: the 12 categorical
# descriptors with their levels (fruit shape and flesh color from the reference assessment),
# and the 5 quantitative profiles cut into N_CLASSES equal-width classes over their range.
# The covering phase runs until the coverage of descriptor classes reaches CORE_COV_MAX, and
# the thickening phase replaces the entries added after CORE_COV_MIN by the rarest accessions.
# The stratified maximum-dissimilarity selection of the same size is kept as a comparison.
# =============================================================================
D   <- as.matrix(base_run$gower)
eqD <- as.matrix(runs$Equal_all$gower)
ids_all <- rownames(main$x)

# ---- class table for ShinyCore ------------------------------------------------
class_table <- function(x, n_classes = N_CLASSES) {
  out <- data.frame(row.names = rownames(x))
  for (tr in CAT) out[[tr]] <- as.character(cat_modes[rownames(x), tr])
  for (tr in CONT) out[[tr]] <- sc_classes(x[[tr]], n_classes)
  out
}
sc_data <- class_table(main$x)
save_table(data.frame(no = rownames(sc_data), sc_data, check.names = FALSE), "core_class_table")
class_def <- bind_rows(lapply(CONT, function(tr) { r <- range(main$x[[tr]], na.rm = TRUE); w <- diff(r) / N_CLASSES
  data.frame(trait = tr, class = seq_len(N_CLASSES), lower = r[1] + (seq_len(N_CLASSES) - 1) * w, upper = r[1] + seq_len(N_CLASSES) * w,
             n_accessions = as.integer(table(factor(sc_classes(main$x[[tr]]), levels = as.character(seq_len(N_CLASSES))))))
}))
save_table(class_def, "core_quantitative_classes")

# ---- ShinyCore selection -------------------------------------------------------
sc <- shinycore(sc_data)
core_ids <- sc$core; N_CORE <- length(core_ids)
stopifnot(!anyDuplicated(core_ids))
cat("ShinyCore: covering phase", sc$N1, "accessions (CV >=", CORE_COV_MAX, "), N0 =", sc$N0, "(CV >", CORE_COV_MIN, "); core size", N_CORE, "\n")
save_table(data.frame(iteration = seq_along(sc$history), no = sc$history, coverage = sc$coverage_history,
                      phase = ifelse(seq_along(sc$history) <= sc$N0, "kept", "replaced in thickening")), "core_covering_history")
save_table(data.frame(no = names(sc$score), rarity_score = as.numeric(sc$score), in_core = names(sc$score) %in% core_ids), "core_rarity_scores")

core <- data.frame(no = core_ids, order = seq_along(core_ids), source = ifelse(seq_along(core_ids) <= sc$N0, "covering", "thickening")) %>%
  left_join(cluster_tab, by = "no") %>% left_join(petal_audit[, c("no", "flower_type")], by = "no")
core <- cbind(core, round(main$x[core$no, CONT, drop = FALSE], 2), cat_modes[core$no, c("fruit_shape", "flesh_color", "leaf_base_shape")])
save_table(core, "core_collection")

# ---- comparison selection: stratified maximum-dissimilarity set of the same size -------------
strata <- table(base_run$cl); quota <- N_CORE * strata / sum(strata); allocation <- pmax(floor(quota), 1)
while (sum(allocation) < N_CORE) { j <- which.max(quota - allocation); allocation[j] <- allocation[j] + 1 }
while (sum(allocation) > N_CORE) { j <- which.max(ifelse(allocation > 1, allocation - quota, -Inf)); allocation[j] <- allocation[j] - 1 }
farthest_subset <- function(D, ids, n) {
  if (n >= length(ids)) return(ids)
  sub <- D[ids, ids, drop = FALSE]
  if (n == 1) return(ids[which.min(rowMeans(sub))])
  pair <- which(sub == max(sub), arr.ind = TRUE)[1, ]; sel <- ids[pair]
  while (length(sel) < n) { left <- setdiff(ids, sel); scores <- apply(D[left, sel, drop = FALSE], 1, min); sel <- c(sel, left[which.max(scores)]) }
  sel
}
maximin_ids <- unlist(lapply(names(strata), function(g) farthest_subset(D, names(base_run$cl)[base_run$cl == as.integer(g)], allocation[[g]])), use.names = FALSE)

# ---- evaluation metrics ------------------------------------------------------------
sc_metrics <- function(sel) {   # ShinyCore report values for an arbitrary set
  sub <- sc_data[sel, , drop = FALSE]; keep <- vapply(sc_data, sc_count, integer(1)) >= 2
  c(class_coverage = 100 * mean(vapply(names(sc_data)[keep], function(j) sc_count(sub[[j]]) / sc_count(sc_data[[j]]), numeric(1))),
    shannon = mean(vapply(names(sc_data)[keep], function(j) sc_shannon(sub[[j]]), numeric(1))),
    rarity = log10(sum(sc$score[sel])))
}
coverage_metrics <- function(sel) {
  cr <- vapply(CONT, function(tr) { v <- main$x[[tr]]; sv <- v[match(sel, rownames(main$x))]; diff(range(sv, na.rm = TRUE)) / diff(range(v, na.rm = TRUE)) }, numeric(1))
  lev_total <- sum(vapply(cat_modes, function(v) length(unique(na.omit(v))), integer(1)))
  lev_sel   <- sum(vapply(cat_modes[sel, , drop = FALSE], function(v) length(unique(na.omit(v))), integer(1)))
  Ds <- D[sel, sel, drop = FALSE]; diag(Ds) <- NA
  c(sc_metrics(sel),
    mean_range_coverage = 100 * mean(cr), modal_level_coverage = 100 * lev_sel / lev_total,
    weighted_nearest_distance = mean(apply(D[, sel, drop = FALSE], 1, min)), equal_nearest_distance = mean(apply(eqD[, sel, drop = FALSE], 1, min)),
    entry_nearest_entry = mean(apply(Ds, 1, min, na.rm = TRUE)), entry_entry = mean(Ds[upper.tri(Ds)]),
    mean_abs_standardized_mean_difference = mean(vapply(CONT, function(tr) abs(mean(main$x[sel, tr], na.rm = TRUE) - mean(main$x[[tr]], na.rm = TRUE)) / sd(main$x[[tr]], na.rm = TRUE), numeric(1))))
}
observed <- coverage_metrics(core_ids); observed_maximin <- coverage_metrics(maximin_ids)
save_table(data.frame(metric = names(observed), ShinyCore = as.numeric(observed), Maximin = as.numeric(observed_maximin),
                      Collection = c(100, sc$SH_all, log10(sum(sc$score)), 100, 100, 0, 0, NA, NA, 0)), "core_method_comparison")
set.seed(SEED + 3)
random_bench <- bind_rows(lapply(seq_len(B_RANDOM), function(b) {
  simple <- sample(ids_all, N_CORE)
  strat  <- unlist(lapply(names(strata), function(g) sample(names(base_run$cl)[base_run$cl == as.integer(g)], allocation[[g]])), use.names = FALSE)
  bind_rows(data.frame(iteration = b, scheme = "Random", as.list(coverage_metrics(simple))),
            data.frame(iteration = b, scheme = "Stratified random", as.list(coverage_metrics(strat))))
}))
save_table(random_bench, "core_random_benchmarks")
bench_summary <- random_bench %>% pivot_longer(-c(iteration, scheme), names_to = "metric", values_to = "value") %>%
  group_by(scheme, metric) %>% summarise(median = median(value), lower = quantile(value, .025), upper = quantile(value, .975), .groups = "drop") %>%
  mutate(core_value = observed[metric], maximin_value = observed_maximin[metric])
save_table(bench_summary, "core_benchmark_summary")
range_tab <- bind_rows(lapply(CONT, function(tr) { v <- main$x[[tr]]; s <- main$x[core_ids, tr]
  data.frame(trait = tr, full_min = min(v, na.rm = TRUE), full_max = max(v, na.rm = TRUE), core_min = min(s, na.rm = TRUE), core_max = max(s, na.rm = TRUE),
             range_coverage_pct = 100 * diff(range(s, na.rm = TRUE)) / diff(range(v, na.rm = TRUE))) }))
save_table(range_tab, "core_continuous_coverage")
level_tab <- bind_rows(lapply(CAT, function(tr) bind_rows(lapply(sort(unique(na.omit(cat_modes[[tr]]))), function(lv)
  data.frame(trait = tr, level = lv, full_n = sum(cat_modes[[tr]] == lv, na.rm = TRUE), core_n = sum(cat_modes[core_ids, tr] == lv, na.rm = TRUE))))))
save_table(level_tab, "core_categorical_coverage")
candidates <- bind_rows(lapply(CONT, function(tr) {
  z <- data.frame(no = rownames(main$x), value = main$x[[tr]]) %>% left_join(years_by_id, by = "no") %>% filter(n_years >= 2)
  z <- if (tr == "bloom_doy") z %>% arrange(value) else z %>% arrange(desc(value))
  z %>% slice_head(n = 5) %>% mutate(trait = tr, in_core = no %in% core_ids)
}))
save_table(candidates, "trait_candidates")

# ---- association of core membership with clusters and observed years --------------------
set.seed(SEED + 4)
chi_core_years <- chisq.test(table(years_by_id$n_years, years_by_id$no %in% core_ids), simulate.p.value = TRUE, B = 10000)
chi_core_cluster <- chisq.test(table(cluster_tab$cluster, cluster_tab$no %in% core_ids), simulate.p.value = TRUE, B = 10000)
core_by_cluster <- data.frame(cluster = as.integer(names(strata)), collection_n = as.integer(strata), core_n = as.integer(table(factor(core[["cluster"]], levels = names(strata)))),
                              maximin_n = as.integer(allocation))
save_table(core_by_cluster, "core_by_cluster")

# ---- sensitivity of the ShinyCore selection to the class definition and the coverage targets -------
sens_core <- bind_rows(lapply(list(list(n = 3L, mx = CORE_COV_MAX, mn = CORE_COV_MIN), list(n = 5L, mx = CORE_COV_MAX, mn = CORE_COV_MIN), list(n = 7L, mx = CORE_COV_MAX, mn = CORE_COV_MIN),
                                   list(n = 5L, mx = 1.00, mn = 0.95), list(n = 5L, mx = 0.95, mn = 0.90)), function(z) {
  r <- shinycore(class_table(main$x, z$n), cov_max = z$mx, cov_min = z$mn); m <- coverage_metrics(r$core)
  data.frame(n_classes = z$n, cov_max = z$mx, cov_min = z$mn, core_size = length(r$core), N0 = r$N0, shared_with_primary = length(intersect(r$core, core_ids)),
             class_coverage = r$CV * 100, shannon = r$SH, mean_range_coverage = m[["mean_range_coverage"]], modal_level_coverage = m[["modal_level_coverage"]],
             entry_nearest_entry = m[["entry_nearest_entry"]])
}))
save_table(sens_core, "core_sensitivity")

write_json(list(n = N_CORE, N0 = sc$N0, N1 = sc$N1, cov_max = CORE_COV_MAX, cov_min = CORE_COV_MIN, n_classes = N_CLASSES, n_markers = sc$n_markers, n_classes_total = sc$n_classes,
                CV = sc$CV, SH = sc$SH, RS = sc$RS, SH_collection = sc$SH_all, RS_collection = log10(sum(sc$score)),
                metrics = as.list(observed), maximin_metrics = as.list(observed_maximin), by_cluster = core_by_cluster, allocation_maximin = as.list(allocation),
                n_levels_total = sum(level_tab$full_n > 0), n_levels_retained = sum(level_tab$core_n > 0), omitted = level_tab[level_tab$core_n == 0, c("trait", "level", "full_n")],
                chi_years = list(statistic = unname(chi_core_years$statistic), p = chi_core_years$p.value),
                chi_cluster = list(statistic = unname(chi_core_cluster$statistic), p = chi_core_cluster$p.value),
                core_by_years = as.list(table(years_by_id$n_years[years_by_id$no %in% core_ids])), n_double = sum(core$flower_type == "double"),
                n_covering = sc$N0, n_thickening = N_CORE - sc$N0, shared_with_maximin = length(intersect(core_ids, maximin_ids)),
                sensitivity = sens_core, coverage_history = sc$coverage_history),
           file.path(OUT, "core_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 5)
