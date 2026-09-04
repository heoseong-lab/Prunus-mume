# =============================================================================
# 02_diversity.R — 계획서 A-1 ~ A-3
#   A-1 기술통계 · Shannon 다양성지수 / A-2 반복성·일치율 / A-3 형질 간 상관
# 선행: 01_preprocess.R
# =============================================================================

if (!exists("acc_mat")) source(file.path("R", "01_preprocess.R"), encoding = "UTF-8")

CAT_TRAITS <- c(trait_of("nominal"), trait_of("ordinal"))
# 자원 매트릭스에서 각 형질에 대응하는 열 이름
acc_col <- function(tr) if (tr %in% trait_of("nominal")) str_c(tr, "__mode") else tr

## ---- A-1 기술통계 -----------------------------------------------------------
dat_long |>
  select(year, all_of(trait_of("continuous")), bloom_doy, dev_days) |>
  pivot_longer(-year, names_to = "trait", values_to = "value") |>
  filter(is.finite(value)) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), CV_pct = 100 * sd / mean,
            min = min(value), median = median(value), max = max(value),
            .by = c(trait, year)) |>
  arrange(trait, year) |>
  save_tab("A1_descriptive_continuous")

dat_long |>
  select(all_of(CAT_TRAITS)) |>
  pivot_longer(everything(), names_to = "trait", values_to = "level") |>
  filter(!is.na(level)) |>
  count(trait, level) |>
  mutate(pct = 100 * n / sum(n), .by = trait) |>
  arrange(trait, desc(n)) |>
  save_tab("A1_descriptive_categorical")

## Shannon–Weaver 다양성지수 (자원 단위, 0~1 정규화)
shannon_tab <- bind_rows(
  tibble(trait = CAT_TRAITS, 유형 = "범주형") |>
    mutate(H_prime  = map_dbl(trait, \(t) shannon(acc_mat[[acc_col(t)]])),
           n_levels = map_int(trait, \(t) n_distinct(acc_mat[[acc_col(t)]], na.rm = TRUE))),
  # 연속형은 모두 연차 보정된 BLUP 을 쓴다 (03 의 군집 분석과 대푯값을 일치시킨다)
  tibble(trait = str_c(c(trait_of("continuous"), "bloom_doy", "dev_days"), "_blup"),
         유형 = "연속형") |>
    filter(trait %in% names(acc_mat)) |>
    mutate(H_prime  = map_dbl(trait, \(t) shannon(acc_mat[[t]])),
           n_levels = NA_integer_)
) |>
  arrange(desc(H_prime)) |>
  save_tab("A1_shannon_diversity")

p_div <- ggplot(shannon_tab, aes(fct_reorder(trait, H_prime), H_prime, fill = 유형)) +
  geom_col(width = .68) +
  scale_fill_manual(
    values = c(범주형 = PAL[["green"]], 연속형 = PAL[["amber"]]),
    name = "Trait type",
    labels = c(범주형 = "Categorical", 연속형 = "Continuous")
  ) +
  scale_x_discrete(labels = trait_label) +
  coord_flip() +
  labs(x = NULL, y = "Shannon diversity index (H')")
save_fig(p_div, "A1_shannon_diversity", w = 8, h = 6.5)

## ---- A-2 반복성 -------------------------------------------------------------
p_rep <- ggplot(repeatability_tab, aes(fct_reorder(trait, repeatability), repeatability, fill = 판정)) +
  geom_col(width = .6) +
  geom_hline(yintercept = c(.2, .5), linetype = "dashed", colour = "grey55") +
  scale_fill_manual(
    values = c(`채택` = PAL[["green"]], `가중치 축소` = PAL[["amber"]],
               `제외 검토` = PAL[["rust"]]),
    name = "Decision",
    labels = c(`채택` = "Retain", `가중치 축소` = "Down-weight",
               `제외 검토` = "Consider exclusion")
  ) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "Repeatability (R)")
save_fig(p_rep, "A2_repeatability", w = 7.5, h = 4)

## 연차 간 상관 — 반복성의 직관적 확인
year_cor <- expand_grid(
    trait = c(trait_of("continuous"), "bloom_doy"),
    pair  = list(c(2024, 2025), c(2025, 2026), c(2024, 2026))
  ) |>
  mutate(res = map2(trait, pair, \(tr, yp) {
    w <- dat_long |>
      filter(year %in% yp) |>
      select(no, year, v = all_of(tr)) |>
      pivot_wider(names_from = year, values_from = v, names_prefix = "y") |>
      drop_na()
    if (nrow(w) < 5) return(tibble(n = nrow(w), r = NA_real_))
    tibble(n = nrow(w), r = cor(w[[2]], w[[3]]))
  })) |>
  unnest(res) |>
  mutate(pair = map_chr(pair, \(p) str_c(p, collapse = "-"))) |>
  save_tab("A2_year_to_year_correlation")

## 범주형 형질의 연차 간 판정 안정성
# 두 지표를 구분한다.
#   concordance_pct  : 모든 조사연차 완전일치 — 형질 자체의 안정성
#   mode_valid_pct   : 과반(2/3 이상) 일치 — '최빈값을 자원 대푯값으로 쓸 수 있는가'
# 자원 대푯값의 신뢰도를 묻는 것이므로 형질 채택 기준은 후자가 맞다.
cat_long <- dat_long |>
  semi_join(acc_mat |> filter(!single_year) |> select(no), by = "no") |>
  select(no, all_of(CAT_TRAITS)) |>
  pivot_longer(-no, names_to = "trait", values_to = "level") |>
  filter(!is.na(level))

concord <- cat_long |>
  summarise(n_obs = n(), max_freq = max(table(level)),
            n_uniq = n_distinct(level), .by = c(trait, no)) |>
  summarise(n_accession = n(),
            concordance_pct = 100 * mean(n_uniq == 1),
            mode_valid_pct  = 100 * mean(max_freq / n_obs > 0.5),
            .by = trait) |>
  left_join(
    cat_long |>
      summarise(kappa = fleiss_kappa(no, level),
                # κ 는 형질별 최다 관측(3개년 완전) 자원만으로 산출된다 — n 을 함께 보고한다.
                n_kappa = fleiss_kappa_n(no, level),
                n_levels = n_distinct(level), .by = trait),
    by = "trait") |>
  # Krippendorff α — 주 지표. 2년 이상 관측 자원 전부(181점)를 쓰고 순서형은 순서를 반영한다.
  left_join(
    map(CAT_TRAITS, \(tr) {
      d   <- filter(cat_long, trait == tr)
      lvl <- if (tr %in% trait_of("ordinal")) "ordinal" else "nominal"
      code <- if (lvl == "ordinal") as.integer(factor(d$level, levels = levels_of(tr)))
              else                  as.integer(factor(d$level))
      tibble(trait = tr, alpha = kripp_alpha(d$no, code, lvl),
             n_alpha = kripp_alpha_n(d$no, code), alpha_metric = lvl)
    }) |> list_rbind(),
    by = "trait") |>
  arrange(desc(alpha)) |>
  mutate(기관 = if_else(str_detect(trait, "^leaf|^flower"), "잎·꽃", "과실·병해")) |>
  save_tab("A2_categorical_concordance")

message("[02] 범주형 연차 일치도 (Krippendorff α 주 지표 / Fleiss κ 대조):")
print(select(concord, trait, n_levels, n_accession, concordance_pct, mode_valid_pct,
             n_alpha, alpha, n_kappa, kappa))

p_conc <- concord |>
  select(trait, organ_group = 기관, complete_agreement = concordance_pct,
         valid_mode_majority = mode_valid_pct, fleiss_kappa = kappa, kripp_alpha = alpha) |>
  mutate(fleiss_kappa = 100 * fleiss_kappa, kripp_alpha = 100 * kripp_alpha) |>
  pivot_longer(-c(trait, organ_group), names_to = "metric", values_to = "pct") |>
  mutate(metric = factor(
    metric,
    levels = c("complete_agreement", "valid_mode_majority", "fleiss_kappa", "kripp_alpha"),
    labels = c("Complete agreement", "Valid mode (majority)",
               "Fleiss' kappa (99 complete accessions)", "Krippendorff's alpha (181 accessions)")
  )) |>
  ggplot(aes(fct_reorder(trait, pct), pct, fill = organ_group)) +
  geom_col(width = .68) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(
    values = c(`잎·꽃` = PAL[["green"]], `과실·병해` = PAL[["rust"]]),
    name = "Organ group",
    labels = c(`잎·꽃` = "Leaf/flower", `과실·병해` = "Fruit/disease")
  ) +
  facet_wrap(~ metric, nrow = 1) +
  coord_flip() +
  labs(x = NULL, y = "Percentage; kappa and alpha multiplied by 100")
save_fig(p_conc, "A2_categorical_concordance", w = 13, h = 5)

## ---- A-2b 순서형 형질의 연차 간 척도 이동 -----------------------------------
# Fleiss κ 는 조사 연차를 교환 가능한 평가자로 본다. 따라서 판정 기준이 해마다
# 통째로 이동하면(예: 2026년에 전반적으로 한 등급 낮게 매김) 그 이동이 '불일치'로
# 나타나 κ 를 끌어내린다. κ 만으로는 "형질에 신호가 없다"와 "잣대가 해마다 달랐다"를
# 구분할 수 없으므로, 두 가지를 따로 계량한다.
#   (1) 3개년 완전조사 99점(균형자료)에서 연차별 평균 UPOV 코드와 그 이동폭
#   (2) 같은 자료에 누적 링크 혼합모형(CLMM)을 적합해 연차 고정효과를 검정
# (1)은 척도 이동의 크기를, (2)는 그 이동이 통계적으로 뒷받침되는지를 답한다.
ORD_TRAITS <- trait_of("ordinal")
ids99_div  <- acc_mat |> filter(n_years == 3) |> pull(no)

year_shift <- ORD_TRAITS |>
  map(\(tr) {
    lv <- TRAITS$levels[[match(tr, TRAITS$trait)]]
    d  <- dat_long |> filter(no %in% ids99_div, !is.na(.data[[tr]])) |>
      mutate(code = match(.data[[tr]], lv))
    if (!nrow(d)) return(NULL)
    m <- summarise(d, mean_code = mean(code), .by = year) |> arrange(year)
    tibble(trait = tr, n_level = length(lv), n_obs = nrow(d),
           !!!set_names(as.list(round(m$mean_code, 2)), str_c("y", m$year)),
           이동폭 = max(m$mean_code) - min(m$mean_code))
  }) |> compact() |> list_rbind()

# CLMM — 있으면 연차 효과를 검정하고, 없으면 이동폭만 보고한다
if (OPT[["ordinal"]]) {
  clmm_tab <- ORD_TRAITS |>
    map(\(tr) {
      lv <- TRAITS$levels[[match(tr, TRAITS$trait)]]
      d  <- dat_long |> filter(no %in% ids99_div) |>
        transmute(accession = factor(no), year_f = factor(year),
                  y = ordered(.data[[tr]], levels = lv)) |>
        filter(!is.na(y)) |> droplevels()
      if (nlevels(d$y) < 2) return(NULL)
      f <- try(ordinal::clmm(y ~ year_f + (1 | accession), data = d, Hess = TRUE),
               silent = TRUE)
      if (inherits(f, "try-error")) return(tibble(trait = tr, ICC = NA_real_, p_year = NA_real_))
      cf <- coef(summary(f)); yr <- cf[grep("^year_f", rownames(cf)), , drop = FALSE]
      vg <- as.numeric(ordinal::VarCorr(f)$accession)
      # 잠재척도 ICC. CLMM 은 잔차 분산을 pi^2/3 으로 고정하므로 자원 간 분리가
      # 뚜렷하면 1 로 포화한다. 형질 서열 확인용이며 Gower 가중치로는 쓰지 않는다.
      tibble(trait = tr, ICC = vg / (vg + pi^2 / 3),
             p_year = if (nrow(yr)) min(yr[, 4]) else NA_real_)
    }) |> compact() |> list_rbind()
  year_shift <- left_join(year_shift, clmm_tab, by = "trait")
} else {
  message("[02] ordinal 미설치 -> 순서형 연차 효과 CLMM 검정 생략 (이동폭만 보고)")
  year_shift <- mutate(year_shift, ICC = NA_real_, p_year = NA_real_)
}

year_shift <- year_shift |>
  left_join(select(concord, trait, kappa, alpha), by = "trait") |>
  arrange(alpha) |>
  save_tab("A2_ordinal_year_shift")

message("[02] 순서형 연차 척도 이동 (3개년 완전조사 ", length(ids99_div), "점 기준):")
print(year_shift)

for (idx in c("alpha", "kappa")) if (sum(is.finite(year_shift[[idx]]) & is.finite(year_shift$이동폭)) >= 5) {
  r_s <- cor(year_shift[[idx]], year_shift$이동폭, method = "spearman", use = "complete.obs")
  neg <- filter(year_shift, .data[[idx]] < 0); pos <- filter(year_shift, .data[[idx]] >= 0)
  message(sprintf(paste0("[02] %s 와 연차 이동폭의 순위상관 = %.3f | ",
                         "%s < 0 형질 평균 이동폭 %.2f등급 (n = %d), >= 0 형질 %.2f등급 (n = %d)"),
                  idx, r_s, idx, mean(neg$이동폭), nrow(neg), mean(pos$이동폭), nrow(pos)))
}
message("     일치도가 낮은 형질일수록 연차 간 판정 척도가 크게 이동했다는 뜻이다.")

## ---- A-3 형질 간 상관 -------------------------------------------------------
# 연속×연속 Spearman | 범주×범주 Cramér's V | 혼합 sqrt(η²)
#
# [세 지표를 한 행렬에 놓는 근거]
# 공통 척도는 "부호 있는 상관"이 아니라 "상관의 크기"이다. Cramér's V 는 2x2 표에서
# |φ| 와 같고, √η² 는 연속형을 범주의 가변수에 회귀했을 때의 다중상관 R(상관비 η)이다.
# 즉 셋 다 [0, 1] 의 상관계수형 크기이다. V 와 √η² 는 방향이 정의되지 않으므로
# (범주에는 고유한 순서가 없다), 척도를 맞추는 방향은 무부호 지표에 부호를 만드는 것이
# 아니라 Spearman 의 부호를 떼는 것이다. 따라서 그림의 색은 |값| 으로 칠하고,
# 부호는 라벨에만 남긴다. 부호는 연속×연속 조합에서만 의미를 갖는다.
#
# [편향보정이 필요한 이유]
# 범위를 맞춰도 영점이 다르다. V 와 η² 는 범주 수가 크고 표본이 작을수록 위로 편향되어
# 독립일 때에도 기댓값이 0 보다 크다(이 자료에서 V 의 귀무 기댓값 중앙값 ≈ 0.12 로,
# 관측 중앙값 0.123 과 사실상 같다). 반면 Spearman 은 귀무 기댓값이 0 이다.
# 그래서 V 는 Bergsma(2013), η² 는 조정 η² = 1 - (1-η²)(n-1)/(n-k) 로 보정한 값을
# 함께 산출하고, 그림과 본문에는 보정값을 쓴다. Spearman 은 보정하지 않는다.
ASSOC_USE_BC <- TRUE     # FALSE 로 두면 그림·행렬이 무보정 원값으로 돌아간다

cont_cols <- str_c(c(trait_of("continuous"), "bloom_doy", "dev_days"), "_blup") |>
  keep(\(x) x %in% names(acc_mat))
cat_cols  <- map_chr(CAT_TRAITS, acc_col) |> keep(\(x) x %in% names(acc_mat))
all_cols  <- c(cont_cols, cat_cols)

# 무보정 / 보정 두 벌을 같은 방식으로 계산한다 (bc = TRUE 면 보정 함수를 쓴다)
assoc_value <- function(a, b, k, bc) {
  if (a == b) return(1)
  switch(k,
    "Spearman"   = suppressWarnings(cor(acc_mat[[a]], acc_mat[[b]],
                                        method = "spearman", use = "pairwise.complete.obs")),
    "Cramér's V" = if (bc) cramers_v_bc(acc_mat[[a]], acc_mat[[b]])
                   else    cramers_v(   acc_mat[[a]], acc_mat[[b]]),
    {
      f <- if (bc) eta_squared_adj else eta_squared
      sqrt(if (a %in% cont_cols) f(acc_mat[[a]], acc_mat[[b]])
           else                  f(acc_mat[[b]], acc_mat[[a]]))
    })
}

assoc_long <- expand_grid(x = all_cols, y = all_cols) |>
  mutate(
    kind     = case_when(x %in% cont_cols & y %in% cont_cols ~ "Spearman",
                         !(x %in% cont_cols) & !(y %in% cont_cols) ~ "Cramér's V",
                         .default = "sqrt(η²)"),
    value    = pmap_dbl(list(x, y, kind), \(a, b, k) assoc_value(a, b, k, bc = FALSE)),
    value_bc = pmap_dbl(list(x, y, kind), \(a, b, k) assoc_value(a, b, k, bc = TRUE)),
    shown    = if (ASSOC_USE_BC) value_bc else value
  )
save_tab(assoc_long, "A3_association_long")

assoc_long |>
  select(x, y, shown) |>
  pivot_wider(names_from = y, values_from = shown) |>
  save_tab("A3_association_matrix")

# 지표별 요약 — 본문 3.4 에 인용하는 수치
assoc_summary <- assoc_long |>
  filter(x != y) |>
  mutate(pair = map2_chr(x, y, \(a, b) str_c(sort(c(a, b)), collapse = "|"))) |>
  distinct(pair, .keep_all = TRUE) |>
  summarise(n              = n(),
            중앙값_무보정  = median(abs(value),    na.rm = TRUE),
            중앙값_보정    = median(abs(value_bc), na.rm = TRUE),
            보정후_0인쌍   = sum(abs(value_bc) < 1e-9, na.rm = TRUE),
            .by = kind) |>
  save_tab("A3_association_summary")
print(assoc_summary)
message("[02] A-3 연관행렬: 색은 |값|, 표시값은 ",
        if (ASSOC_USE_BC) "편향보정 후" else "무보정", " 기준")

# 색은 크기 |값| 으로 (무부호 지표가 92% 이므로 발산형 척도는 오독을 부른다),
# 라벨은 부호를 살려 둔다 — 부호는 연속×연속 조합에만 있다.
p_heat <- ggplot(assoc_long, aes(x, y, fill = abs(shown))) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = if_else(abs(shown) >= .35 & x != y, sprintf("%.2f", shown), "")),
            size = 2.4, colour = "grey15") +
  scale_fill_gradient(low = "grey96", high = PAL[["green"]],
                      limits = c(0, 1), na.value = "grey92",
                      name = "Association strength") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7.5), panel.grid = element_blank())
save_fig(p_heat, "A3_association_heatmap", w = 10, h = 9)

message("[02] A-1 ~ A-3 완료")
