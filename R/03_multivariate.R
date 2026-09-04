# =============================================================================
# 03_multivariate.R — 계획서 A-4 ~ A-7
#   A-4 형질 신뢰도 · FAMD · 가중 Gower PCoA
#   A-5 군집 (주 분석 + 민감도 6종 + ARI 비교)
#   A-6 핵심집단   A-7 목표별 선발
# 선행: 02_diversity.R
#
# 설계 원칙 — 형질을 버리지 않는다.
#   20형질을 모두 넣되 Gower 거리에서 형질별 신뢰도로 가중한다.
#   신뢰도가 0 인 형질은 자동으로 기여가 0 이 되므로, 임의의 컷오프를 두는 대신
#   자료가 스스로 가중치를 정하게 한다.
# =============================================================================

if (!exists("concord")) source(file.path("R", "02_diversity.R"), encoding = "UTF-8")

## ---- A-4a 형질 신뢰도 -------------------------------------------------------
# 주 분석 가중치 (2026-09-02 결정 C):
#   연속형: 혼합모형 반복성 R (연차·배치 보정) → Spearman-Brown(R, n̄) = 자원 BLUP 의 신뢰도
#   범주형: Fleiss κ (연차를 평가자로 본 우연 보정 일치도, 3개년 완전조사 자원)
# 두 값은 서로 다른 구성개념이다(연속형은 n̄ 회 대푯값의 신뢰도, 범주형은 단일 평정 간
# 일치도). 둘 다 [0, 1] 의 재현성 계수로서 '상대 가중치'로만 쓰며 동치라고 주장하지 않는다.
#
# 대조 지표 두 가지를 같은 표에 싣는다:
#   (a) BLUP신뢰도_PEV — 자원별 1 - PEV_i/σ²G 의 평균. S-B(n̄) 근사가 실제 예측값의 신뢰도와
#       얼마나 가까운지 보인다(이 자료에서 최대 차 0.04).
#   (b) alpha / alpha_SB — Krippendorff α(순서형은 순서 반영, 2년 이상 181점)와 그 S-B 환산값.
#       α 를 주 가중치로 쓰지 않는 이유: α + 대칭 S-B 를 적용하면 명목형 엽저 형태의 가중치가
#       0.458 → 0.780 으로 올라, 명목형의 0/1 부분거리(1.0)와 곱해져 단일 형질 기여 0.78 이
#       겹꽃 대비(0.5 × 0.99 = 0.495)를 넘어선다. 그 결과 주 군집이 엽저 형태 하나로 결정된다
#       (구성 ⑧ 에서 그대로 보여 준다). 이는 재현성의 발견이 아니라 Gower 의 명목형 코딩이
#       만드는 지렛대 효과이므로, 주 분석은 κ 를 유지하고 α 는 민감도 구성으로 제시한다.
N_BAR <- mean(acc_mat$n_years)

reliability <- bind_rows(
  repeatability_tab |>
    transmute(trait, 유형 = "연속형", 지표 = "Spearman-Brown(R, n̄)",
              원지표 = repeatability, 신뢰도 = spearman_brown(repeatability, N_BAR),
              BLUP신뢰도_PEV = blup_reliability_pev, alpha = NA_real_, alpha_SB = NA_real_),
  concord |>
    transmute(trait, 유형 = "범주형", 지표 = "Fleiss κ",
              원지표 = kappa, 신뢰도 = kappa,
              BLUP신뢰도_PEV = NA_real_, alpha = alpha, alpha_SB = spearman_brown(alpha, N_BAR))
) |>
  mutate(가중치 = pmax(신뢰도, 0),
         등급 = case_when(신뢰도 >= .60 ~ "높음",
                          신뢰도 >= .40 ~ "보통",
                          신뢰도 >= .20 ~ "낮음",
                          .default = "우연 수준")) |>
  arrange(desc(신뢰도)) |>
  save_tab("A4_trait_reliability")

message("[03] 형질 신뢰도 (연속형 Spearman-Brown n̄=", round(N_BAR, 2), " / 범주형 Fleiss κ; 대조: PEV, α):")
print(select(reliability, trait, 유형, 원지표, 신뢰도, 가중치, 등급, BLUP신뢰도_PEV, alpha, alpha_SB), n = 25)
message("[03] 연속형: S-B(n̄) 와 PEV 기반 BLUP 신뢰도의 최대 차 = ",
        round(max(abs(reliability$신뢰도 - reliability$BLUP신뢰도_PEV), na.rm = TRUE), 3))

p_rel <- ggplot(reliability, aes(fct_reorder(trait, 신뢰도), 신뢰도, fill = 등급)) +
  geom_col(width = .68) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(
    values = c(`높음` = PAL[["green"]], `보통` = PAL[["sage"]],
               `낮음` = PAL[["amber"]], `우연 수준` = PAL[["rust"]]),
    name = "Reliability class",
    labels = c(`높음` = "High", `보통` = "Moderate", `낮음` = "Low",
               `우연 수준` = "Chance-level")
  ) +
  coord_flip() +
  facet_wrap(~ 유형, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(연속형 = "Continuous", 범주형 = "Categorical"))) +
  labs(x = NULL, y = "Reliability (weight)")
save_fig(p_rel, "A4_trait_reliability", w = 8, h = 7)

## ---- 분석용 프레임 ----------------------------------------------------------
mv_cont <- reliability |> filter(유형 == "연속형") |> pull(trait) |>
  str_c("_blup") |> keep(\(x) x %in% names(acc_mat))
mv_nom  <- trait_of("nominal") |> str_c("__mode") |> keep(\(x) x %in% names(acc_mat))
mv_ord  <- trait_of("ordinal")  |> str_c("_code")  |> keep(\(x) x %in% names(acc_mat))
mv_all  <- c(mv_cont, mv_nom, mv_ord)

# 열 이름 -> 형질 이름 -> 가중치
col_to_trait <- \(x) str_remove(x, "(_blup|__mode|_code)$")
w_all <- reliability$가중치[match(col_to_trait(mv_all), reliability$trait)] |> replace_na(0)
names(w_all) <- mv_all
# 대안 가중 체계(연속형 S-B(R), 범주형 S-B(α)) — 민감도 구성 ⑧ 용 (리뷰 권고안)
w_alpha <- reliability |>
  mutate(w = if_else(유형 == "연속형", pmax(신뢰도, 0), pmax(alpha_SB, 0))) |>
  pull(w, name = trait)
w_alpha <- w_alpha[col_to_trait(mv_all)] |> replace_na(0); names(w_alpha) <- mv_all
message("[03] 구성 ⑧ 가중치에서 leaf_base_shape = ", round(w_alpha[["leaf_base_shape__mode"]], 3),
        " (주 분석 ", round(w_all[["leaf_base_shape__mode"]], 3), ")")

mv_x <- acc_mat |>
  select(no, all_of(mv_all)) |>
  mutate(across(all_of(mv_nom), \(x) factor(replace_na(x, "unknown")))) |>
  column_to_rownames("no")

message("[03] 분석 매트릭스 ", nrow(mv_x), " x ", ncol(mv_x), "형질 (전 형질 포함)")

# 어느 형질이 거리를 지배하는가
w_share <- tibble(형질 = col_to_trait(mv_all), 가중치 = as.numeric(w_all)) |>
  mutate(기여율_pct = round(100 * 가중치 / sum(가중치), 1)) |>
  arrange(desc(가중치)) |>
  save_tab("A4_weight_share")
message("[03] 상위 3형질 누적 기여율 ",
        sum(head(w_share$기여율_pct, 3)), "% | 가중치 0 형질 ",
        sum(w_share$가중치 == 0), "개")

# 20형질 완전 관측 자원 (사용자 제안 검증용)
complete_ids <- dat_long |>
  mutate(ok = if_all(all_of(c(trait_of("continuous"), trait_of("nominal"),
                              trait_of("ordinal"), "bloom_doy", "maturity_doy")),
                     \(x) !is.na(x))) |>
  summarise(any_complete = any(ok), .by = no) |>
  filter(any_complete) |> pull(no)
message("[03] 20형질 완전 관측 자원: ", length(complete_ids), " / ", nrow(mv_x))

# 3개년 모두 조사된 자원. 연차 간 재현성 추정의 근거 집단이자,
# "완전 균형 자료만 써야 한다"는 관행적 요구를 그대로 따랐을 때 남는 표본이다.
# 이 요구를 따르면 105점(2개년 82 + 1개년 23)을 버리게 되므로, 그 대가로 무엇을
# 얻는지 자원 순위와 군집 구조 양쪽에서 확인한다.
ids99 <- acc_mat |> filter(n_years == 3) |> pull(no)
message("[03] 3개년 완전조사 자원: ", length(ids99), " / ", nrow(mv_x),
        " (제외되는 자원 ", nrow(mv_x) - length(ids99), "점)")

## ---- A-4b FAMD (전 20형질, 형질 기여도 확인용) ------------------------------
famd <- FAMD(mv_x, ncp = 5, graph = FALSE)

eig <- as_tibble(famd$eig, rownames = "dim") |>
  set_names("dim", "eigenvalue", "pct_var", "cum_pct_var") |>
  save_tab("A4_famd_eigenvalues")

famd_ind <- as_tibble(famd$ind$coord, rownames = "no") |>
  left_join(select(acc_mat, no, n_years), by = "no") |>
  save_tab("A4_famd_individual_coords")

save_fig(fviz_screeplot(famd, addlabels = TRUE, ncp = 8,
                        barfill = PAL[["green"]], barcolor = PAL[["green"]]) +
           labs(x = "Dimension", y = "Explained variance (%)"),
         "A4_famd_scree", w = 7, h = 4.5)

save_fig(fviz_famd_var(famd, repel = TRUE, col.var = PAL[["green"]]) +
           labs(title = NULL),
         "A4_famd_variables", w = 7.5, h = 6.5)

## ---- A-4c 가중 Gower 거리 (주 분석) -----------------------------------------
# [거리의 기하 — 왜 제곱근을 취하는가]
# daisy(metric = "gower") 는 1 - s 를 돌려준다. Gower(1971)가 증명한 것은 s 가 양반정치라서
# sqrt(1 - s) 가 유클리드 공간에 내장된다는 것이다. 즉 1 - s 는 유클리드 거리의 '제곱'이다.
# Ward 기준은 유클리드 거리 위에서 정의되고, hclust(method = "ward.D2") 는 입력을 다시
# 제곱해서 쓰므로(Murtagh & Legendre 2014), 1 - s 를 그대로 넣으면 유클리드 거리의 4제곱에
# Ward 를 적용하는 셈이 된다. cmdscale(PCoA) 도 같은 이유로 1 - s 에서는 음의 고유값이
# 대량으로 나온다(이 자료에서 204개 중 153개). sqrt(1 - s) 로 바꾸면 음의 고유값이 7개로
# 줄고(합 비율 0.27%), Ward.D2·PCoA·실루엣이 모두 같은 유클리드 거리 위에서 돌아간다.
# 이 변환은 단조이므로 핵심집단의 maximin 선발(A-6)에는 영향이 없다.
#
# 검증(2026-09-02): ward.D on (1-s) 와 ward.D2 on sqrt(1-s) 는 k = 2..4 에서 ARI = 1.000 으로
# 동치. 이전 구현(ward.D2 on 1-s)과는 k = 2 에서 ARI -0.07 로 첫 분할이 달랐다.
gower_sq <- function(x, weights = NULL) {
  d <- if (is.null(weights)) daisy(x, metric = "gower")
       else                  daisy(x, metric = "gower", weights = weights)
  sqrt(d)
}
gower_w <- gower_sq(mv_x, w_all)                        # 신뢰도 가중, 유클리드형
# 균등 가중 Gower 는 민감도 구성(runs$E 등)에서 make_dist() 로 다시 계산한다.

# 유클리드성 진단 — 변환 전후의 음의 고유값을 표로 남긴다 (방법 2.8 의 근거)
eig_check <- function(d, label) {
  e <- cmdscale(d, k = attr(d, "Size") - 1, eig = TRUE)$eig
  tibble(거리 = label, n_eig = length(e),
         n_negative = sum(e < -1e-10), pct_negative = round(100 * mean(e < -1e-10), 1),
         neg_over_pos_pct = round(100 * abs(sum(e[e < 0])) / sum(e[e > 0]), 2),
         min_eig = round(min(e), 4), max_eig = round(max(e), 4))
}
bind_rows(eig_check(daisy(mv_x, metric = "gower", weights = w_all), "1 - s (raw Gower)"),
          eig_check(gower_w, "sqrt(1 - s)")) |>
  save_tab("A4_pcoa_eigen_check") |> print()

pcoa <- cmdscale(gower_w, k = 5, eig = TRUE)
# sqrt(1 - s) 는 사실상 유클리드라 음의 고유값이 거의 없다. 분모는 관례대로 양의 합.
pcoa_var <- 100 * pcoa$eig[1:5] / sum(pcoa$eig[pcoa$eig > 0])

pcoa_ind <- tibble(no = rownames(mv_x), PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2]) |>
  left_join(select(acc_mat, no, n_years), by = "no")

## 교락 점검 — 주축이 유전변이가 아니라 수확 배치를 재현하고 있지 않은가
# 자원마다 연차별 배치가 각 1회씩만 나타나 "최빈값"이 성립하지 않는다.
# 실제로는 최초 조사 연차의 배치가 되므로 이름에 그대로 밝힌다.
batch_first <- dat_long |>
  filter(!is.na(harvest_batch)) |>
  arrange(no, year) |>
  summarise(batch_first = as.character(first(harvest_batch)),
            maturity_doy = mean(maturity_doy, na.rm = TRUE), .by = no)

confound <- pcoa_ind |>
  left_join(batch_first, by = "no") |>
  left_join(select(acc_mat, no, bloom_doy_blup), by = "no")

# eta2 는 군 수가 많으면 우연만으로도 커진다. 귀무 기대값 (g-1)/(n-1) 을 함께 싣는다.
c("PCo1", "PCo2") |>
  map(\(d) {
    ok <- !is.na(confound$batch_first) & is.finite(confound[[d]])
    g  <- n_distinct(confound$batch_first[ok])
    tibble(축 = d,
           eta2_harvest_batch = eta_squared(confound[[d]], confound$batch_first),
           eta2_null_expected = (g - 1) / (sum(ok) - 1),
           n_batch = g,
           r_maturity_doy = suppressWarnings(cor(confound[[d]], confound$maturity_doy,
                             method = "spearman", use = "pairwise.complete.obs")),
           r_bloom_doy = suppressWarnings(cor(confound[[d]], confound$bloom_doy_blup,
                             method = "spearman", use = "pairwise.complete.obs")))
  }) |>
  list_rbind() |> save_tab("A4_confounding_check") |> print()

## ---- A-5 군집: 주 분석 + 민감도 5종 -----------------------------------------
# 거리 정의와 연결법을 함께 바꾸어 비교한다.
#   gower_w  신뢰도 가중 Gower (주 분석)
#   gower    균등 가중 Gower
#   euclid_z 전 형질을 수치로 바꾼 뒤 z 표준화 + 유클리드
#            — 유전자원 표현형 군집의 관행적 절차(Marklová et al. 2026 등)를 재현한 것.
#              표준화는 형질별 분산을 1로 맞추므로 재현성과 무관하게 균등 가중이 된다.
best_k <- function(d, link = "ward.D2", ks = 2:10) {
  h <- hclust(d, method = link)
  s <- map_dbl(ks, \(k) mean(silhouette(cutree(h, k), d)[, 3]))
  lst(hc = h, k = ks[which.max(s)], sil = max(s),
      scan = tibble(k = ks, silhouette = s))
}

make_dist <- function(x, traits, kind) {
  # Gower 계열은 모두 sqrt(1 - s) (A-4c 참조). 유클리드 계열은 그대로.
  switch(kind,
    gower_w  = gower_sq(x, w_all[traits]),
    gower_wa = gower_sq(x, w_alpha[traits]),         # 대안 가중(S-B(R) + S-B(α))
    gower    = gower_sq(x),
    euclid_z = {
      # 주의 — 명목형을 as.numeric() 으로 바꾸면 수준의 나열 순서(알파벳 순)가
      # 그대로 거리에 들어간다. 임의의 선택이므로 코딩 민감도를 별도로 계량한다
      # (A5_euclid_coding_sensitivity).
      z <- as.data.frame(map(x, \(v) {
             v <- as.numeric(v)                       # 명목형은 수준 코드로
             replace(v, is.na(v), mean(v, na.rm = TRUE)) }))
      rownames(z) <- rownames(x)
      dist(scale(z), method = "euclidean")
    },
    stop("알 수 없는 거리 종류: ", kind))
}

runs <- list(
  W = lst(label = "\u2460 Reliability-weighted 20 traits (primary analysis)", ids = rownames(mv_x), traits = mv_all,
          dist_kind = "gower_w",  link = "ward.D2"),
  E = lst(label = "\u2461 Equal-weighted 20 traits", ids = rownames(mv_x), traits = mv_all,
          dist_kind = "gower",    link = "ward.D2"),
  C = lst(label = "\u2462 Complete-observation accessions, 20 traits", ids = complete_ids, traits = mv_all,
          dist_kind = "gower",    link = "ward.D2"),
  S = lst(label = "\u2463 Composite reliability \u2265 .4 (Spearman-Brown S or \u03ba)", ids = rownames(mv_x),
          traits = mv_all[w_all >= .4], dist_kind = "gower", link = "ward.D2"),
  Z = lst(label = "\u2464 Standardized Euclidean distance, 20 traits", ids = rownames(mv_x), traits = mv_all,
          dist_kind = "euclid_z", link = "ward.D2"),
  U = lst(label = "\u2465 Equal-weighted 20 traits + UPGMA", ids = rownames(mv_x), traits = mv_all,
          dist_kind = "gower",    link = "average"),
  T3 = lst(label = "\u2466 Accessions surveyed in all three years, reliability-weighted",
           ids = ids99, traits = mv_all, dist_kind = "gower_w", link = "ward.D2"),
  A  = lst(label = "\u2467 Krippendorff \u03b1 with Spearman-Brown on both data types",
           ids = rownames(mv_x), traits = mv_all, dist_kind = "gower_wa", link = "ward.D2")
)

# 축 눈금용 짧은 라벨. 그림에서 제목·부제를 없앴으므로 구성 기호의 뜻을
# 축 자체가 설명해야 한다 (가운뎃점은 글꼴 의존이라 쉼표로 쓴다).
RUN_SHORT <- c(
  W = "\u2460 Reliability-weighted\nGower, Ward.D2",
  E = "\u2461 Equal-weighted\nGower, Ward.D2",
  C = "\u2462 Complete cases\nGower, Ward.D2",
  S = "\u2463 Reliability >= 0.4\nGower, Ward.D2",
  Z = "\u2464 Standardized\nEuclidean, Ward.D2",
  U  = "\u2465 Equal-weighted\nGower, UPGMA",
  T3 = "\u2466 Three-year complete\nWeighted Gower, Ward.D2",
  A  = "\u2467 Krippendorff \u03b1\n+ S-B on both types")

runs <- imap(runs, \(r, nm) {
  x <- mv_x[r$ids, r$traits, drop = FALSE] |> droplevels()
  d <- make_dist(x, r$traits, r$dist_kind)
  b <- best_k(d, r$link)
  c(r, lst(n = nrow(x), n_trait = ncol(x), dist = d, hc = b$hc, k = b$k,
           sil = b$sil, scan = b$scan,
           cl = tibble(no = rownames(x), cluster = cutree(b$hc, b$k))))
})

run_summary <- imap(runs, \(r, nm) tibble(구성 = nm, 설명 = r$label, n자원 = r$n,
                                          n형질 = r$n_trait, 거리 = r$dist_kind, 연결 = r$link,
                                          k = r$k, silhouette = round(r$sil, 3),
                                          군집크기 = str_c(sort(as.integer(table(r$cl$cluster)),
                                                               decreasing = TRUE),
                                                           collapse = " / "))) |>
  list_rbind() |> save_tab("A5_run_comparison")
message("[03] 군집 구성 비교:")
print(run_summary)

## 구성 ⑦ 보강 — 99점만으로 BLUP 을 다시 추정해도 자원 순위가 유지되는가
# 군집 구조(ARI)만으로는 "자원을 버려도 되는가"에 답할 수 없다. 선발의 근거가 되는
# 자원 순위가 유지되는지도 보아야 하므로, 99점만으로 같은 혼합모형을 다시 적합하여
# 전체 204점 BLUP 과의 순위상관을 계산한다.
blup99_rank <- c(trait_of("continuous"), "bloom_doy") |>
  map(\(tr) {
    d <- dat_long |> filter(no %in% ids99) |>
      select(accession, year_f, harvest_batch, y = all_of(tr)) |> filter(is.finite(y))
    if (n_distinct(d$accession) < 10) return(NULL)
    form <- if (tr %in% trait_of("continuous") && n_distinct(d$harvest_batch) > 2) {
      y ~ year_f + (1 | accession) + (1 | harvest_batch)
    } else {
      y ~ year_f + (1 | accession)
    }
    f <- try(lmer(form, data = d,
                  control = lmerControl(check.conv.singular = "ignore")), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    b <- ranef(f)$accession |> rownames_to_column("no") |> as_tibble() |>
      transmute(no, blup_99 = `(Intercept)`)
    m <- acc_mat |>
      select(no, blup_all = all_of(str_c(tr, "_blup"))) |>
      filter(no %in% ids99) |>
      inner_join(b, by = "no") |>
      filter(is.finite(blup_all))
    tibble(형질 = tr, n = nrow(m),
           spearman = cor(m$blup_99, m$blup_all, method = "spearman"),
           pearson  = cor(m$blup_99, m$blup_all))
  }) |>
  compact() |> list_rbind() |>
  save_tab("A5_subset99_blup_rank")
message("[03] 99점 재추정 BLUP vs 전체 204점 BLUP 순위상관 = ",
        str_c(sprintf("%.3f", range(blup99_rank$spearman)), collapse = " ~ "))
print(blup99_rank)

## ARI — 구성 간 군집 구조가 얼마나 같은가
RUN_ORD <- names(runs)
ari_tab <- expand_grid(a = RUN_ORD, b = RUN_ORD) |>
  mutate(ARI = map2_dbl(a, b, \(i, j) {
    m <- inner_join(runs[[i]]$cl, runs[[j]]$cl, by = "no", suffix = c("_a", "_b"))
    if (nrow(m) < 10) return(NA_real_)
    adj_rand(m$cluster_a, m$cluster_b)
  })) |>
  save_tab("A5_ARI_between_runs")
message("[03] 구성 간 ARI (1 = 동일 구조, 0 = 우연 수준):")
print(pivot_wider(ari_tab, names_from = b, values_from = ARI))

## 구성 ⑤ 의 코딩 민감도 — 명목형을 정수로 바꾸는 순간 임의 선택이 들어간다
# 표준화 유클리드 절차는 명목형 형질을 수치로 바꾸어야 하는데, 그 정수 코드는
# 수준의 나열 순서(관행적으로 알파벳 순)에 따라 정해질 뿐 아무 의미가 없다.
# 순서를 바꾸면 결과가 얼마나 달라지는지 직접 계량한다.
N_PERM <- 200
set.seed(2026)
euclid_perm <- map(seq_len(N_PERM), \(i) {
  x <- mv_x
  for (v in mv_nom) x[[v]] <- factor(x[[v]], levels = sample(levels(x[[v]])))
  d <- make_dist(x, mv_all, "euclid_z")
  b <- best_k(d, "ward.D2")
  m <- inner_join(runs$W$cl, tibble(no = rownames(x), cluster = cutree(b$hc, b$k)),
                  by = "no", suffix = c("_a", "_b"))
  tibble(perm = i, k = b$k, silhouette = b$sil,
         ARI_vs_primary = adj_rand(m$cluster_a, m$cluster_b))
}) |> list_rbind() |> save_tab("A5_euclid_coding_sensitivity")

euclid_sens <- tibble(
  n_perm = N_PERM,
  ARI_observed = ari_tab |> filter(a == "W", b == "Z") |> pull(ARI),
  ARI_min = min(euclid_perm$ARI_vs_primary), ARI_max = max(euclid_perm$ARI_vs_primary),
  ARI_median = median(euclid_perm$ARI_vs_primary), ARI_sd = sd(euclid_perm$ARI_vs_primary),
  k_2 = sum(euclid_perm$k == 2), k_3 = sum(euclid_perm$k == 3),
  k_4plus = sum(euclid_perm$k >= 4)) |>
  save_tab("A5_euclid_coding_sensitivity_summary")
print(euclid_sens)
message(sprintf(paste0("[03] 구성 5 는 명목형 코딩 순서에 좌우된다: ARI(1,5) = %.3f ~ %.3f ",
                       "(중앙값 %.3f, SD %.3f), k = 2/3/4+ 가 %d/%d/%d 회"),
                euclid_sens$ARI_min, euclid_sens$ARI_max, euclid_sens$ARI_median,
                euclid_sens$ARI_sd, euclid_sens$k_2, euclid_sens$k_3, euclid_sens$k_4plus))

## 형질 하나의 효과 — 과육색 6수준의 순서 720가지 전수 열거 (나머지 3형질은 알파벳 순 고정)
# 방법 2.9 에 기술된 계산. 재현성이 우연 수준(κ < 0)인 형질 하나의 표기 순서만으로
# 구성 ⑤ 의 결론이 얼마나 흔들리는지 보인다.
flesh_col <- "flesh_color__mode"
flesh_lv  <- levels(mv_x[[flesh_col]])
perm_all  <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), \(i) lapply(perm_all(v[-i]), \(p) c(v[i], p))))
flesh_perm <- imap(perm_all(flesh_lv), \(ord, i) {
  x <- mv_x; x[[flesh_col]] <- factor(x[[flesh_col]], levels = ord)
  b <- best_k(make_dist(x, mv_all, "euclid_z"), "ward.D2")
  m <- inner_join(runs$W$cl, tibble(no = rownames(x), cluster = cutree(b$hc, b$k)),
                  by = "no", suffix = c("_a", "_b"))
  tibble(perm = i, order = str_c(ord, collapse = ">"), k = b$k, silhouette = b$sil,
         ARI_vs_primary = adj_rand(m$cluster_a, m$cluster_b))
}) |> list_rbind() |> save_tab("A5_euclid_flesh_color_enumeration")
flesh_sens <- tibble(n_orderings = nrow(flesh_perm),
                     ARI_min = min(flesh_perm$ARI_vs_primary), ARI_max = max(flesh_perm$ARI_vs_primary),
                     ARI_median = median(flesh_perm$ARI_vs_primary),
                     k_min = min(flesh_perm$k), k_max = max(flesh_perm$k)) |>
  save_tab("A5_euclid_flesh_color_summary")
message(sprintf("[03] 과육색 순서 720가지: ARI(1,5) = %.3f ~ %.3f (중앙값 %.3f), k = %d ~ %d",
                flesh_sens$ARI_min, flesh_sens$ARI_max, flesh_sens$ARI_median, flesh_sens$k_min, flesh_sens$k_max))

lab_of <- \(n) str_c(n, " = ", runs[[n]]$label)
p_ari <- ari_tab |>
  mutate(a = factor(a, levels = RUN_ORD), b = factor(b, levels = rev(RUN_ORD))) |>
  ggplot(aes(a, b, fill = ARI)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("%.2f", ARI)), size = 3.4, colour = "grey12") +
  scale_fill_gradient2(low = PAL[["rust"]], mid = "grey96", high = PAL[["green"]],
                       midpoint = 0.3, limits = c(-.1, 1), name = "Adjusted Rand index") +
  scale_x_discrete(labels = RUN_SHORT) +
  scale_y_discrete(labels = RUN_SHORT) +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5, lineheight = .95),
        axis.text.y = element_text(size = 7.5, lineheight = .95))
save_fig(p_ari, "A5_ARI_heatmap", w = 10.5, h = 8.2)

## 주 분석 채택
PRIMARY <- runs$W
K_BEST  <- PRIMARY$k
cluster_tab <- PRIMARY$cl |> rename(cluster_ward = cluster) |> save_tab("A5_cluster_assignment")
acc_mat <- left_join(acc_mat, cluster_tab, by = "no")

save_fig(ggplot(PRIMARY$scan, aes(k, silhouette)) +
           geom_line(colour = "grey60") +
           geom_point(aes(colour = k == K_BEST), size = 3) +
           scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = PAL[["rust"]]), guide = "none") +
           scale_x_continuous(breaks = 2:10) +
           labs(x = "Number of clusters (k)", y = "Mean silhouette width"),
         "A5_silhouette_scan", w = 6.5, h = 4)

p_pcoa <- pcoa_ind |>
  left_join(cluster_tab, by = "no") |>
  mutate(cluster_ward = factor(cluster_ward)) |>
  ggplot(aes(PCo1, PCo2)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  stat_ellipse(aes(colour = cluster_ward), linewidth = .4, alpha = .6) +
  geom_point(aes(colour = cluster_ward), size = 2, alpha = .85) +
  # 라벨 배치를 재현 가능하게 만드는 두 가지. 둘 다 있어야 한다.
  #  seed     : ggrepel 은 기본적으로 그리는 시점의 난수 상태를 쓴다.
  #  max.time : 기본값 0.5초로, 반복 계산을 "시간"으로 끊는다. 기계 부하에 따라
  #             완료되는 반복 횟수가 달라져 같은 자료로도 라벨 위치가 바뀐다.
  #             Inf 로 두면 max.iter(기본 10000) 만으로 끊겨 결정적이 된다.
  geom_text_repel(aes(label = no), size = 1.9, max.overlaps = 10,
                  seed = 2026, max.time = Inf, max.iter = 10000,
                  segment.colour = "grey78", colour = "grey35") +
  labs(x = sprintf("PCoA1 (%.1f%%)", pcoa_var[1]), y = sprintf("PCoA2 (%.1f%%)", pcoa_var[2]),
       colour = "Cluster")
save_fig(p_pcoa, "A4_pcoa_biplot", w = 10, h = 8.5)

save_fig(fviz_dend(PRIMARY$hc, k = K_BEST, cex = .38, horiz = TRUE, rect = TRUE,
                   rect_border = "grey60", k_colors = scales::hue_pal()(K_BEST)) +
           labs(title = NULL), 
         "A5_dendrogram", w = 9, h = 15)

## 군집별 프로파일 — 전 20형질 (신뢰도 병기)
acc_mat |>
  select(cluster_ward, all_of(mv_cont)) |>
  pivot_longer(-cluster_ward, names_to = "col", values_to = "v") |>
  summarise(n = sum(is.finite(v)), mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
            .by = c(cluster_ward, col)) |>
  mutate(trait = col_to_trait(col),
         신뢰도 = round(reliability$신뢰도[match(trait, reliability$trait)], 3)) |>
  save_tab("A5_cluster_profile_continuous")

acc_mat |>
  select(cluster_ward, all_of(mv_nom)) |>
  pivot_longer(-cluster_ward, names_to = "col", values_to = "level") |>
  count(cluster_ward, col, level) |>
  mutate(pct = 100 * n / sum(n), .by = c(cluster_ward, col)) |>
  slice_max(n, n = 2, by = c(cluster_ward, col), with_ties = FALSE) |>
  mutate(trait = col_to_trait(col),
         신뢰도 = round(reliability$신뢰도[match(trait, reliability$trait)], 3)) |>
  save_tab("A5_cluster_profile_categorical")

## 주 분석 군집을 무엇이 가르는가 — 형질별 판별력
cluster_driver <- bind_rows(
  tibble(형질 = col_to_trait(mv_cont), col = mv_cont) |>
    mutate(지표 = "η²",
           값 = map_dbl(col, \(v) eta_squared(acc_mat[[v]], acc_mat$cluster_ward))),
  tibble(형질 = col_to_trait(c(mv_nom, mv_ord)), col = c(mv_nom, mv_ord)) |>
    mutate(지표 = "Cramér's V",
           값 = map_dbl(col, \(v) cramers_v(acc_mat[[v]], acc_mat$cluster_ward)))
) |>
  left_join(select(reliability, 형질 = trait, 신뢰도), by = "형질") |>
  arrange(desc(값)) |>
  select(형질, 지표, 판별력 = 값, 신뢰도) |>
  save_tab("A5_cluster_drivers")

message("[03] 주 분석 군집 판별력 상위:")
print(head(cluster_driver, 6))

## ---- A-6 핵심집단 -----------------------------------------------------------
CORE_FRAC <- 0.15
gm <- as.matrix(gower_w)   # sqrt(1-s). 단조 변환이라 maximin 선발 결과는 1-s 와 같다.

pick_maxmin <- function(ids, n_pick) {
  if (n_pick >= length(ids)) return(ids)
  sub <- gm[ids, ids, drop = FALSE]
  # n_pick = 1 이면 최대거리 쌍에서 시작할 수 없다 — 가장 바깥에 있는 1점을 고른다.
  if (n_pick == 1L) return(ids[which.max(rowSums(sub))])
  ij  <- which(sub == max(sub), arr.ind = TRUE)[1, ]
  sel <- ids[c(ij[[1]], ij[[2]])]
  while (length(sel) < n_pick) {
    rest <- setdiff(ids, sel)
    sel  <- c(sel, rest[which.max(apply(gm[rest, sel, drop = FALSE], 1, min))])
  }
  sel
}

core_ids <- acc_mat |>
  filter(no %in% rownames(mv_x)) |>
  summarise(core = list(pick_maxmin(no, max(1, round(n() * CORE_FRAC)))), .by = cluster_ward) |>
  pull(core) |> list_c()

# 대표성 지표의 정의를 문헌 표준에 맞춘다.
#   MD% / VD% : 평균·분산이 전체집단과 유의하게 다른 형질의 비율 (수치형)
#   CR%       : Hu et al.(2000) 의 coincidence rate of range — 수치형 형질의 범위 유지율.
#               범주형에는 정의되지 않으므로 범주 수준 유지율을 따로 계산해 병기한다.
#   VR%       : variable rate of coefficient of variation — 변이계수의 비.
#               (이전 구현은 분산비였다. 참고용으로 분산비도 함께 싣는다.)
num_repr <- function(cols, cd) {
  map(cols, \(v) {
    a <- keep(acc_mat[[v]], is.finite); b <- keep(cd[[v]], is.finite)
    if (length(b) < 3 || length(unique(a)) < 2) return(NULL)
    rng_a <- diff(range(a))
    tibble(trait = col_to_trait(v), col = v,
           n_core = length(b), n_total = length(a),
           md = t.test(b, a)$p.value < .05,
           vd = var.test(b, a)$p.value < .05,
           cv_ratio  = (sd(b) / mean(b)) / (sd(a) / mean(a)),
           var_ratio = var(b) / var(a),
           range_cr  = if (rng_a > 0) diff(range(b)) / rng_a else NA_real_)
  }) |> list_rbind()
}
lvl_cr <- function(cols, cd) map_dbl(cols, \(v) {
  lw <- unique(na.omit(acc_mat[[v]])); lc <- unique(na.omit(cd[[v]]))
  length(intersect(lw, lc)) / length(lw)
})

core_eval <- {
  cd <- filter(acc_mat, no %in% core_ids)
  cs <- num_repr(mv_cont, cd)                 # 연속형 6형질
  os <- num_repr(mv_ord,  cd)                 # 순서형 10형질 (UPOV 코드)
  cr_nom <- lvl_cr(mv_nom, cd)                # 명목형 4형질
  cr_ord <- lvl_cr(mv_ord, cd)                # 순서형 10형질

  bind_rows(mutate(cs, 형질군 = "연속형"), mutate(os, 형질군 = "순서형")) |>
    save_tab("A6_core_representativeness_by_trait")

  tibble(n_core = length(core_ids), n_total = nrow(mv_x),
         비율_pct = 100 * length(core_ids) / nrow(mv_x),
         MD_pct = 100 * mean(cs$md), VD_pct = 100 * mean(cs$vd),
         MD_all_pct = 100 * mean(c(cs$md, os$md)),
         VD_all_pct = 100 * mean(c(cs$vd, os$vd)),
         CR_range_pct = 100 * mean(cs$range_cr, na.rm = TRUE),
         CR_level_pct = 100 * mean(c(cr_nom, cr_ord)),
         CR_level_nominal_pct = 100 * mean(cr_nom),
         VR_pct = 100 * mean(cs$cv_ratio),
         VR_variance_pct = 100 * mean(cs$var_ratio))
} |> save_tab("A6_core_representativeness")
print(core_eval)

acc_mat |>
  filter(no %in% core_ids) |>
  select(no, cluster_ward, n_years, all_of(c(mv_cont, mv_nom))) |>
  arrange(cluster_ward, no) |>
  save_tab("A6_core_collection")

## ---- A-7 육종 목표별 선발 ---------------------------------------------------
top_by <- function(col, n = 15, desc = TRUE, label) {
  if (!col %in% names(acc_mat)) return(NULL)
  d <- filter(acc_mat, is.finite(.data[[col]]))
  # 동률로 목록이 늘어나지 않도록 자원번호를 2차 정렬키로 두고 정확히 n 점만 뽑는다.
  d <- arrange(d, no)
  d <- if (desc) slice_max(d, .data[[col]], n = n, with_ties = FALSE)
       else      slice_min(d, .data[[col]], n = n, with_ties = FALSE)
  rel <- reliability$신뢰도[match(col_to_trait(col), reliability$trait)]
  d |> transmute(목표 = label, no, cluster = cluster_ward, 형질 = col,
                 값 = round(.data[[col]], 3), n_years, 신뢰도 = round(rel, 3),
                 판정 = case_when(is.na(rel) ~ "미평가", rel >= .6 ~ "높음",
                                  rel >= .4 ~ "보통", rel >= .2 ~ "낮음 — 재확인 권장",
                                  .default = "우연 수준 — 근거로 쓰지 말 것"))
}

list(
  top_by("fruit_weight_blup",           15, TRUE,  "대과성"),
  top_by("soluble_solids_content_blup", 15, TRUE,  "고당도"),
  top_by("titratable_acidity_blup",     15, TRUE,  "고산 (가공적성)"),
  top_by("titratable_acidity_blup",     15, FALSE, "저산 (생식용)"),
  top_by("petal_number_blup",           15, TRUE,  "다판성 (관상)"),
  top_by("bloom_doy_blup",              15, FALSE, "조생"),
  top_by("bloom_doy_blup",              15, TRUE,  "만생")
) |> compact() |> list_rbind() |> save_tab("A7_selection_by_objective")

dat_long |>
  summarise(n_obs = n(),
            n_br = sum(!is.na(brown_rot_severity)), n_sh = sum(!is.na(shot_hole_severity)),
            # all(logical(0)) 은 TRUE 다 — 관측이 실제로 있는 자원만 무병으로 본다.
            brown_rot_free = n_br >= 2 && all(brown_rot_severity == "0%", na.rm = TRUE),
            shot_hole_free = n_sh >= 2 && all(shot_hole_severity == "0%", na.rm = TRUE),
            .by = no) |>
  filter(n_obs >= 2, brown_rot_free, shot_hole_free) |>
  save_tab("A7_disease_free_accessions")

message("[03] A-4 ~ A-7 완료 | 주 분석 = 가중 20형질, k = ", K_BEST)
