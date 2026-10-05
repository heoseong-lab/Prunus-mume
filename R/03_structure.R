# =============================================================================
# 03_structure.R  Weighted Gower distance, ordination, clustering, sensitivity
# =============================================================================
euclideanize <- function(D) {          # Lingoes correction when sqrt(Gower) is not Euclidean
  D <- as.matrix(D); n <- nrow(D); H <- diag(n) - matrix(1 / n, n, n)
  ev <- eigen(-.5 * H %*% (D^2) %*% H, symmetric = TRUE, only.values = TRUE)$values
  cadd <- max(0, -min(ev)); if (cadd < 1e-10) cadd <- 0
  Dc <- sqrt(pmax(D^2 + 2 * cadd, 0)); diag(Dc) <- 0; dimnames(Dc) <- dimnames(D)
  list(distance = as.dist(Dc), constant = cadd, negative_fraction = abs(sum(ev[ev < -1e-10])) / sum(ev[ev > 0]))
}
make_distance <- function(x, w) {
  cols <- intersect(names(w), names(x)); w <- w[cols]; keep <- is.finite(w) & w > 0
  stopifnot(any(keep)); x <- x[, cols[keep], drop = FALSE]; w <- w[keep]
  d <- suppressWarnings(daisy(x, metric = "gower", weights = unname(w)))
  if (any(!is.finite(d))) stop("No comparable positive-weight components in at least one accession pair")
  out <- euclideanize(sqrt(as.matrix(d))); out$gower <- d; out
}
cluster_run <- function(x, w, k = NULL, link = "ward.D2") {
  dis <- make_distance(x, w); h <- hclust(dis$distance, method = link); ks <- 2:min(10, nrow(x) - 1)
  sil <- vapply(ks, function(j) mean(silhouette(cutree(h, j), dis$distance)[, 3]), numeric(1))
  best <- ks[which.max(sil)]; if (is.null(k)) k <- best
  cl <- cutree(h, k); pm <- tapply(x$petal_number, cl, mean, na.rm = TRUE); ord <- order(pm, decreasing = TRUE)
  map <- setNames(seq_along(ord), names(pm)[ord]); cl <- setNames(as.integer(map[as.character(cl)]), names(cl))
  list(cl = cl, k = k, best_k = best, sil = mean(silhouette(cl, dis$distance)[, 3]), scan = data.frame(k = ks, silhouette = sil),
       hc = h, d = dis$distance, gower = dis$gower, correction = dis$constant, negative_fraction = dis$negative_fraction)
}

primary_w <- main$weights[PRIMARY]
base_run  <- cluster_run(main$x, primary_w); K <- base_run$k
cat("Primary k =", K, " silhouette =", round(base_run$sil, 4), " sizes =", paste(table(base_run$cl), collapse = "/"), "\n")
save_table(base_run$scan, "silhouette_scan")
cluster_tab <- data.frame(no = names(base_run$cl), cluster = unname(base_run$cl)) %>% left_join(years_by_id, by = "no")
save_table(cluster_tab, "cluster_assignments")

# ---- sensitivity configurations ------------------------------------------------
equal_w <- setNames(as.numeric(!(PRIMARY %in% FIXED_TRAITS)), PRIMARY)   # equal weights for the measured variables
rel_by  <- function(col) setNames(main$reliability[[col]], main$reliability$trait)
configs <- list(
  Primary            = list(x = main$x, w = primary_w),
  Equal_all          = list(x = main$x, w = equal_w),
  Equal_retained     = list(x = main$x, w = setNames(as.numeric(primary_w > 0), PRIMARY)),
  Alpha_categorical  = list(x = main$x, w = replace(primary_w, setdiff(CAT, FIXED_TRAITS), pmax(0, rel_by("alpha")[setdiff(CAT, FIXED_TRAITS)]))),
  Complete_profiles  = list(x = main$x[complete.cases(main$x[, PRIMARY]), ], w = primary_w))
for (nyr in c(2, 3)) {
  ids <- years_by_id$no[years_by_id$n_years >= nyr]; a <- profile_and_weights(dat[dat$no %in% ids, ])
  configs[[paste0("Min_", nyr, "_years")]] <- list(x = a$x, w = a$weights[PRIMARY])
}

runs <- lapply(configs, function(z) cluster_run(z$x, z$w))
sens <- bind_rows(lapply(names(runs), function(nm) {
  z <- runs[[nm]]; shared <- intersect(names(base_run$cl), names(z$cl)); fixed <- cutree(z$hc, K)
  data.frame(configuration = nm, n_accessions = length(z$cl), n_positive_weights = sum(configs[[nm]]$w > 0), best_k = z$k,
             silhouette = z$sil, cluster_sizes = paste(as.integer(table(z$cl)), collapse = "/"),
             ARI_selected_k = ari(base_run$cl[shared], z$cl[shared]), ARI_fixed_k = ari(base_run$cl[shared], fixed[shared]),
             Lingoes_constant = z$correction, negative_eigen_fraction = z$negative_fraction)
}))
save_table(sens, "sensitivity_configurations")
save_table(bind_rows(lapply(names(configs), function(nm) data.frame(configuration = nm, trait = names(configs[[nm]]$w), weight = unname(configs[[nm]]$w)))), "weights_by_configuration")

# ---- descriptive profiles of the clusters ------------------------------------
profile_c <- bind_rows(lapply(CONT, function(tr) data.frame(no = rownames(main$x), value = main$x[[tr]]) %>%
  left_join(cluster_tab, by = "no") %>% group_by(cluster) %>%
  summarise(n = sum(is.finite(value)), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), .groups = "drop") %>% mutate(trait = tr)))
save_table(profile_c, "cluster_continuous_profiles")
cat_modes <- as.data.frame(setNames(lapply(CAT, function(tr) vapply(rownames(main$x), function(i) mode_value(dat[[tr]][dat$no == i]), character(1))), CAT))
rownames(cat_modes) <- rownames(main$x)
profile_cat <- bind_rows(lapply(CAT, function(tr) data.frame(no = rownames(cat_modes), level = cat_modes[[tr]]) %>%
  left_join(cluster_tab, by = "no") %>% count(cluster, level) %>% group_by(cluster) %>% mutate(percent = 100 * n / sum(n), trait = tr) %>% ungroup()))
save_table(profile_cat, "cluster_categorical_profiles")

# flower type (derived descriptor) and its cross-classification with leaf base shape and cluster
petal_audit <- dat %>% group_by(no) %>% summarise(n_obs = sum(!is.na(petal_number)), minimum = min(petal_number), mean = mean(petal_number),
                                                  maximum = max(petal_number), .groups = "drop") %>%
  mutate(flower_type = ifelse(mean >= DOUBLE_MIN, "double", "single")) %>% left_join(cluster_tab, by = "no")
save_table(petal_audit, "petal_number_accession_audit")
xtab <- petal_audit %>% mutate(leaf_base = cat_modes[no, "leaf_base_shape"]) %>% count(cluster, flower_type, leaf_base)
save_table(xtab, "cluster_by_flower_type_leaf_base")
two_way <- petal_audit %>% mutate(leaf_base = cat_modes[no, "leaf_base_shape"]) %>%
  mutate(group = paste(flower_type, leaf_base, sep = "_"))
ari_two_way <- ari(two_way$cluster, two_way$group)
# agreement of the partition with each categorical descriptor and with selected two-way classifications
ari_cat <- bind_rows(lapply(CAT, function(tr) { ok <- !is.na(cat_modes[[tr]])
  data.frame(classification = tr, ARI = ari(cluster_tab$cluster[match(rownames(cat_modes)[ok], cluster_tab$no)], cat_modes[[tr]][ok])) }))
ari_cat <- bind_rows(ari_cat,
  data.frame(classification = "flower_type", ARI = ari(petal_audit$cluster, petal_audit$flower_type)),
  data.frame(classification = "flower_type x leaf_base_shape", ARI = ari_two_way),
  data.frame(classification = "fruit_shape x flesh_color", ARI = ari(cluster_tab$cluster[match(rownames(cat_modes), cluster_tab$no)], paste(cat_modes$fruit_shape, cat_modes$flesh_color))),
  data.frame(classification = "flower_type x fruit_shape", ARI = ari(petal_audit$cluster, paste(petal_audit$flower_type, cat_modes[petal_audit$no, "fruit_shape"]))))
save_table(ari_cat, "cluster_agreement_with_descriptors")
# the partition obtained with the reference-assessed descriptors at full weight, against fruit shape x flesh color
xtab_fruit <- petal_audit %>% mutate(fruit_shape = cat_modes[no, "fruit_shape"], flesh_color = cat_modes[no, "flesh_color"]) %>% count(cluster, fruit_shape, flesh_color)
save_table(xtab_fruit, "cluster_by_fruit_shape_flesh_color")
set.seed(SEED)
ny <- table(cluster = cluster_tab$cluster, n_years = cluster_tab$n_years)
chi_years <- chisq.test(ny, simulate.p.value = TRUE, B = 10000)
save_table(as.data.frame.matrix(ny) %>% mutate(cluster = rownames(.)), "cluster_by_n_years")

# ---- weight shares versus empirical distance contributions -------------------
x <- main$x[, PRIMARY]; n <- nrow(x); den <- matrix(0, n, n); pieces <- list()
for (tr in PRIMARY) {
  v <- x[[tr]]; ok <- outer(!is.na(v), !is.na(v), "&")
  if (is.factor(v)) dj <- outer(as.character(v), as.character(v), "!=") * 1 else {
    ran <- diff(range(v, na.rm = TRUE)); dj <- if (ran > 0) abs(outer(v, v, "-")) / ran else matrix(0, n, n) }
  dj[!ok] <- 0; pieces[[tr]] <- primary_w[[tr]] * dj; den <- den + primary_w[[tr]] * ok
}
contrib <- vapply(pieces, function(m) sum((m / pmax(den, 1e-15))[upper.tri(m)]), numeric(1))
weight_share <- data.frame(trait = PRIMARY, weight = primary_w, weight_share_pct = 100 * primary_w / sum(primary_w), contribution_pct = 100 * contrib / sum(contrib))
save_table(weight_share, "weight_and_distance_contributions")

# ---- ordination ---------------------------------------------------------------
pc <- cmdscale(base_run$d, k = 5, eig = TRUE); pc_pct <- 100 * pc$eig[1:5] / sum(pc$eig[pc$eig > 0])
pc_tab <- data.frame(no = rownames(pc$points), PCo1 = pc$points[, 1], PCo2 = pc$points[, 2]) %>% left_join(cluster_tab, by = "no") %>%
  left_join(petal_audit[, c("no", "flower_type")], by = "no")
save_table(pc_tab, "pcoa_coordinates")
save_table(data.frame(axis = seq_along(pc$eig), eigenvalue = pc$eig), "pcoa_eigenvalues")
stopifnot(length(PRIMARY) == 17, all(is.finite(base_run$d)), abs(sum(weight_share$contribution_pct) - 100) < 1e-8)
