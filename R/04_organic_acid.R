# =============================================================================
# 04_organic_acid.R — 계획서 S-1 ~ S-8 : 유기산 자원 간 비교통계 (89자원)
#
# 전제: 반복 3회는 서로 다른 개체(생물학적 반복). 오차항이 생물학적 변이이므로
#       자원 간 검정은 유전자원 수준의 추론으로 성립한다.
#
# 입력 우선순위
#   1) 유기산_반복원자료.csv : no, rep, oxalic_acid, malic_acid, citric_acid
#   2) 매실_유기산_자원.csv  : mean/sd 로부터 반복치 재구성
#
#   재구성: n=3, 표본 SD 가 s 일 때 (m-s, m, m+s) 는 평균 m 과 표본 SD s 를
#   정확히 재현한다. ANOVA·LSD·Tukey·Bartlett 은 군평균·군분산·n 만의 함수이므로
#   결과가 원자료와 일치한다. 반면 순위 기반 검정과 잔차 정규성 진단은 무의미하므로
#   원자료가 있을 때만 실행한다.
# =============================================================================

if (!exists("PATH")) source(file.path("R", "00_setup.R"), encoding = "UTF-8")

ACIDS     <- c("oxalic_acid", "malic_acid", "citric_acid")
ACID_MEAN <- c("oxalic_mean", "malic_mean", "citric_mean")
DERIVED   <- c("total_acid", "citric_malic_ratio")
N_REP     <- 3L
ALPHA     <- 0.05

## ---- 보고 단위 --------------------------------------------------------------
# 원자료(매실_유기산_자원.csv)는 mg·(100 g)^-1 FW 로 기록되어 있다.
# 매실 유기산 문헌의 관행(Kang et al. 2020; Yu et al. 2015)에 맞추어
# 보고 단위를 mg·g^-1 FW 로 통일한다.
#
# 자료를 읽는 시점에 농도 열을 한 번만 환산하면 이후 통계량은 자동으로 맞는다.
#   평균·SD·분위수·최소·최대·LSD·HSD·SE·신뢰구간 : x OA_SCALE     (1차 동차)
#   평균제곱 MS                                   : x OA_SCALE^2   (2차 동차)
#   F·p·Bartlett·Welch·CV%·RSD%·왜도·첨도·상관·
#   구성비(%)·구연산/사과산 비·Shapiro p·TOST p   : 척도 불변 — 손대지 않는다
# (Shapiro-Wilk 은 위치·척도 불변이고, log 변환은 상수 log(100) 만큼의
#  평행이동이므로 정규성 진단 결과도 바뀌지 않는다.)
OA_SCALE    <- 1 / 100
OA_UNIT_TXT <- "mg/g FW"          # 콘솔 메시지·표 주석용 (ASCII)

# 축·범례 라벨은 plotmath 로 쓴다. 유니코드 위첨자(⁻¹)는 글꼴에 따라 깨지지만
# plotmath 는 그래픽 장치가 직접 조판하므로 어느 환경에서나 같게 나온다.
OA_LAB        <- expression("Concentration (mg" %.% g^-1 * " FW)")
OA_LAB_OXALIC <- expression("Oxalic acid (mg" %.% g^-1 * " FW)")
OA_LAB_MALIC  <- expression("Malic acid (mg" %.% g^-1 * " FW)")
OA_LAB_CITRIC <- expression("Citric acid (mg" %.% g^-1 * " FW)")

# 패싯 제목(strip)용 — label_parsed 로 조판하므로 plotmath 문법의 문자열이다.
ACID_EXPR <- c(
  oxalic_mean        = "'Oxalic acid (mg'%.%g^-1*' FW)'",
  malic_mean         = "'Malic acid (mg'%.%g^-1*' FW)'",
  citric_mean        = "'Citric acid (mg'%.%g^-1*' FW)'",
  total_acid         = "'Total acid (mg'%.%g^-1*' FW)'",
  citric_malic_ratio = "'Citric/malic ratio'",
  oxalic_acid        = "'Oxalic acid'",
  malic_acid         = "'Malic acid'",
  citric_acid        = "'Citric acid'")
relabel_acid <- function(x) factor(x, levels = names(ACID_EXPR), labels = ACID_EXPR)

# 환산 대상은 농도 열뿐이다. *_rsd_pct · *_pct_of_total · citric_malic_ratio 는
# 비(比)라서 불변이고, ta_mean·ssc_mean·fruit_weight_mean 은 다른 형질이다.
CONC_COLS <- c("oxalic_mean", "oxalic_sd", "malic_mean", "malic_sd",
               "citric_mean", "citric_sd", "total_acid")

## ---- S-1 자료 확정 ----------------------------------------------------------
oa <- read_csv(PATH$oa, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) |>
  mutate(across(any_of(CONC_COLS), \(x) x * OA_SCALE))
message("[04] 농도 단위 환산 mg/(100 g) -> ", OA_UNIT_TXT,
        " | 환산 열 ", sum(CONC_COLS %in% names(oa)), " / ", length(CONC_COLS))

RAW_AVAILABLE <- file.exists(PATH$oa_raw)

rep_dat <- if (RAW_AVAILABLE) {
  message("[04] 반복 원자료 사용: ", basename(PATH$oa_raw))
  read_csv(PATH$oa_raw, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) |>
    mutate(no = as.character(no), rep = as.integer(rep),
           across(any_of(ACIDS), \(x) x * OA_SCALE))   # 원자료도 같은 단위로
} else {
  message("[04] 반복 원자료 없음 -> mean±sd 에서 n=3 반복치 재구성")
  oa |>
    select(no, matches("^(oxalic|malic|citric)_(mean|sd)$")) |>
    pivot_longer(-no, names_to = c("acid", ".value"),
                 names_pattern = "^(oxalic|malic|citric)_(mean|sd)$") |>
    rename(m = mean, s = sd) |>
    uncount(N_REP, .id = "rep") |>
    mutate(value = m + s * c(-1, 0, 1)[rep], acid = str_c(acid, "_acid")) |>
    select(no, rep, acid, value) |>
    pivot_wider(names_from = acid, values_from = value)
}

# QC 결정
EXCLUDE_CITRIC_OUTLIER <- TRUE     # KJM182 citric 550.1 mg/g (차순위의 8.4배)
outlier_ids <- oa |> filter(str_detect(qc_flag, "citric_outlier")) |> pull(no)
highrsd_ids <- oa |> filter(str_detect(qc_flag, "high_rsd")) |> pull(no)
keep_ids    <- if (EXCLUDE_CITRIC_OUTLIER) setdiff(oa$no, outlier_ids) else oa$no

oa_use  <- filter(oa, no %in% keep_ids)
rep_use <- rep_dat |> filter(no %in% keep_ids) |> mutate(no = factor(no))

message("[04] S-1 대상 ", nrow(oa_use), "자원 | 제외 ", str_c(outlier_ids, collapse = ", "),
        " | high_rsd ", length(intersect(highrsd_ids, keep_ids)), "점 포함")

## ---- S-2 기술통계와 분포 ----------------------------------------------------
long_mean <- oa_use |>
  select(no, all_of(c(ACID_MEAN, DERIVED))) |>
  pivot_longer(-no, names_to = "trait", values_to = "value")

desc_oa <- long_mean |>
  summarise(n = n(), mean = mean(value), sd = sd(value), CV_pct = 100 * sd / mean,
            min = min(value), Q1 = quantile(value, .25), median = median(value),
            Q3 = quantile(value, .75), max = max(value),
            skew = skewness(value), kurt = kurtosis(value),
            shapiro_p = shapiro.test(value)$p.value,
            shapiro_p_log = shapiro.test(log(value))$p.value,
            .by = trait) |>
  mutate(권장변환 = case_when(shapiro_p >= .05 ~ "원자료",
                              shapiro_p_log >= .05 ~ "로그 변환",
                              .default = "비모수 검정")) |>
  save_tab("S2_organic_acid_descriptive")

print(select(desc_oa, trait, n, mean, sd, CV_pct, skew, shapiro_p, 권장변환))

oa_use |>
  select(no, ends_with("_rsd_pct")) |>
  pivot_longer(-no, names_to = "trait", values_to = "rsd") |>
  summarise(median_rsd = median(rsd), max_rsd = max(rsd), n_over20 = sum(rsd > 20),
            .by = trait) |>
  save_tab("S2_replicate_rsd")

# 패싯마다 단위가 다르다(농도 4개 + 무차원 비 1개). 공통 x축 라벨로는
# 표현할 수 없으므로 단위를 패싯 제목에 싣는다.
long_lab <- mutate(long_mean, trait = relabel_acid(trait))

p_hist <- ggplot(long_lab, aes(value)) +
  geom_histogram(bins = 22, fill = PAL[["green"]], colour = "white", linewidth = .25) +
  facet_wrap(~ trait, scales = "free", ncol = 3, labeller = label_parsed) +
  labs(x = NULL, y = "Number of accessions")
save_fig(p_hist, "S2_distribution", w = 10, h = 6)

p_qq <- ggplot(long_lab, aes(sample = value)) +
  stat_qq(size = .9, colour = PAL[["green"]]) +
  stat_qq_line(colour = PAL[["rust"]]) +
  facet_wrap(~ trait, scales = "free", ncol = 3, labeller = label_parsed) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles")
save_fig(p_qq, "S2_qqplot", w = 10, h = 6)

## ---- S-3 전체 유의성 검정 ---------------------------------------------------
# H0: mu_1 = ... = mu_k  (자원 간 함량 차이 없음)
anova_tab <- ACIDS |>
  map(\(a) {
    d   <- rep_use |> select(no, y = all_of(a)) |> filter(is.finite(y))
    fit <- aov(y ~ no, data = d)
    s   <- summary(fit)[[1]]
    dfe <- s[["Df"]][2]; mse <- s[["Mean Sq"]][2]
    k   <- n_distinct(d$no)
    tibble(acid = a,
           df_between = s[["Df"]][1], df_within = dfe,
           MS_between = s[["Mean Sq"]][1], MS_within = mse,
           F_value = s[["F value"]][1], p_value = s[["Pr(>F)"]][1],
           bartlett_p = bartlett.test(y ~ no, data = d)$p.value,
           welch_p = oneway.test(y ~ no, data = d, var.equal = FALSE)$p.value,
           LSD_05 = qt(1 - ALPHA/2, dfe) * sqrt(2 * mse / N_REP),
           HSD_05 = qtukey(1 - ALPHA, k, dfe) * sqrt(mse / N_REP))
  }) |>
  list_rbind() |>
  save_tab("S3_anova_summary")

print(anova_tab)

if (RAW_AVAILABLE) {
  ACIDS |>
    map(\(a) {
      kt <- kruskal.test(reformulate("no", a), data = rep_use)
      tibble(acid = a, KW_chisq = kt$statistic, KW_df = kt$parameter, KW_p = kt$p.value)
    }) |> list_rbind() |> save_tab("S3_kruskal_wallis")
} else {
  message("[04] 재구성 자료 -> Kruskal-Wallis / 잔차 정규성 검정 생략 (원자료 필요)")
}

## ---- S-4 평균 분리 : 3갈래 --------------------------------------------------
## (a) 함량 등급화 — 전 자원이 한 장에 들어가는 5등급표
GRADES <- c("매우 낮음", "낮음", "중간", "높음", "매우 높음")
oa_use |>
  select(no, all_of(c(ACID_MEAN, DERIVED))) |>
  mutate(across(-no, \(x) cut(x, quantile(x, seq(0, 1, .2)),
                              include.lowest = TRUE, labels = GRADES),
                .names = "{.col}_grade")) |>
  save_tab("S4a_content_grades")

# 등급 폭과 LSD 의 관계는 두 방향을 구분해야 한다.
#   폭 >= LSD : 한 등급 이상 떨어진 자원은 통계적으로 구분된다 (등급이 실험 해상도보다 거칠다)
#   폭 <= LSD : 같은 등급 안의 자원은 통계적으로 구분되지 않는다
# 두 명제는 동치가 아니다. 최소 등급 폭이 LSD 를 넘으면 성립하는 것은 전자이며,
# 같은 등급 안에서도 최대 (등급 폭) 만큼 차이가 날 수 있으므로 후자는 성립하지 않는다.
grade_check <- ACIDS |>
  map(\(a) {
    v   <- oa_use[[str_replace(a, "_acid$", "_mean")]]
    lsd <- anova_tab |> filter(acid == a) |> pull(LSD_05)
    qs  <- quantile(v, seq(0, 1, .2))
    wid <- diff(qs)
    tibble(acid = a, LSD = lsd,
           min_grade_width = min(wid), max_grade_width = max(wid),
           등급간_구분가능 = min(wid) >= lsd,     # 인접 등급끼리도 통계적으로 구분된다
           등급내_구분불가 = max(wid) <= lsd,     # 같은 등급이면 구분되지 않는다
           해석 = if_else(min(wid) >= lsd,
                          "등급이 다르면 통계적으로 구분됨 (등급 해상도 > LSD)",
                          "일부 인접 등급은 통계적으로 구분되지 않음"))
  }) |> list_rbind() |> save_tab("S4a_grade_width_vs_LSD")
print(grade_check)
message("[04] S-4a 최소 등급폭 >= LSD -> '등급이 다르면 구분됨'이 성립. ",
        "'같은 등급 = 구분 불가'는 폭 <= LSD 일 때만 성립하며 본 자료에서는 성립하지 않는다.")

## (b) 표적 비교 — 상위 10 · 하위 10 자원만 Tukey HSD
tukey_tab <- ACIDS |>
  map(\(a) {
    mcol <- str_replace(a, "_acid$", "_mean")
    ord  <- arrange(oa_use, desc(.data[[mcol]]))
    ids  <- c(head(ord$no, 10), tail(ord$no, 10))
    d    <- rep_use |> filter(no %in% ids) |> mutate(no = droplevels(no)) |>
      select(no, y = all_of(a))
    # 오차항은 부분집합이 아니라 전체 88자원 ANOVA 의 풀링 MSE(df = 176)를 쓴다.
    # S-3·S-4c 와 오차항을 일치시키고 검정력도 유지된다.
    row  <- filter(anova_tab, acid == a)
    hsd  <- agricolae::HSD.test(y = d$y, trt = d$no,
                                DFerror = row$df_within, MSerror = row$MS_within,
                                group = TRUE)
    # HSD.test 를 벡터 인자로 호출하면 평균 열 이름이 호출식에서 만들어지므로
    # 이름이 아니라 위치로 집는다 (1열 = 평균, 2열 = 문자군).
    g <- hsd$groups |> rownames_to_column("no") |> as_tibble()
    names(g)[2:3] <- c("mean_y", "group_y")
    g |>
      transmute(acid = a, 구분 = if_else(no %in% head(ord$no, 10), "상위 10", "하위 10"),
                no, mean = mean_y, group = group_y)
  }) |> list_rbind() |> save_tab("S4b_top_bottom10_tukey")

## (c) 신뢰구간 forest plot — 풀링 MSE 기반 95% CI
ci_tab <- ACIDS |>
  map(\(a) {
    row <- filter(anova_tab, acid == a)
    se  <- sqrt(row$MS_within / N_REP)
    tcr <- qt(1 - ALPHA/2, row$df_within)
    oa_use |>
      transmute(acid = a, no, mean = .data[[str_replace(a, "_acid$", "_mean")]],
                se = se, lwr = mean - tcr * se, upr = mean + tcr * se,
                rank = rank(-mean, ties.method = "first"),
                qc = if_else(no %in% highrsd_ids, "high_rsd", "정상"))
  }) |> list_rbind() |> save_tab("S4c_accession_means_CI")

# 자원 순서를 총산 기준 하나로 통일해 y축 라벨을 한 벌만 쓰고, 폭을 줄여 라벨을 키운다.
acc_order <- oa_use |> arrange(total_acid) |> pull(no)
p_forest <- ci_tab |>
  mutate(no = factor(no, levels = acc_order),
         acid = relabel_acid(fct_relevel(acid, "oxalic_acid", "malic_acid", "citric_acid"))) |>
  ggplot(aes(no, mean)) +
  geom_linerange(aes(ymin = lwr, ymax = upr), colour = "grey55", linewidth = .4) +
  geom_point(aes(colour = qc), size = 1.3) +
  scale_colour_manual(
    values = c(정상 = PAL[["green"]], high_rsd = PAL[["rust"]]),
    name = "QC status",
    labels = c(정상 = "Normal", high_rsd = "High RSD")
  ) +
  facet_wrap(~ acid, scales = "free_x", nrow = 1, labeller = label_parsed) +
  coord_flip() +
  labs(x = NULL, y = OA_LAB) +
  theme(axis.text.y = element_text(size = 5.6),
        panel.grid.major.y = element_line(linewidth = .15, colour = "grey93"),
        legend.position = "top")
save_fig(p_forest, "S4c_forest_plot", w = 9, h = 13)

## ---- S-5 산 조성 유형 -------------------------------------------------------
# 분류 기준은 log(구연산/사과산 비)의 3분위이므로 라벨도 상대적이어야 한다.
# 이 집단은 전체적으로 구연산 우점이라(비 중앙값 1.9) "사과산 우세형"이라는
# 절대 표현은 사실과 어긋난다 — 하위 3분위 30자원 중 21자원이 구연산 > 사과산이다.
# 절대 우점 여부는 별도 열(절대우점)로 함께 보고한다.
CM_LAB <- c("C/M비 하위형", "C/M비 중위형", "C/M비 상위형")
comp <- oa_use |>
  select(no, oxalic = oxalic_pct_of_total, malic = malic_pct_of_total,
         citric = citric_pct_of_total, citric_malic_ratio, total_acid) |>
  mutate(log_cm = log(citric_malic_ratio),
         조성유형 = cut(log_cm,
                        c(-Inf, quantile(log_cm, 1/3), quantile(log_cm, 2/3), Inf),
                        labels = CM_LAB),
         절대우점 = case_when(citric_malic_ratio > 1 ~ "구연산 우점",
                              citric_malic_ratio < 1 ~ "사과산 우점",
                              .default = "동률")) |>
  save_tab("S5_composition_type")
print(count(comp, 조성유형, 절대우점))
message("[04] S-5 라벨은 C/M비 3분위(상대 구분)다. 절대 기준으로 사과산이 우점하는 자원은 ",
        sum(comp$절대우점 == "사과산 우점"), " / ", nrow(comp), "점뿐이다.")

type_dat <- comp |> left_join(select(oa_use, no, all_of(ACID_MEAN)), by = "no")

c(ACID_MEAN, "total_acid") |>
  map(\(v) {
    fit <- aov(reformulate("조성유형", v), data = type_dat)
    s   <- summary(fit)[[1]]
    agricolae::HSD.test(fit, "조성유형", group = TRUE)$groups |>
      rownames_to_column("조성유형") |>
      as_tibble() |>
      transmute(변수 = v, 조성유형, mean = .data[[v]], group = groups,
                F = s[["F value"]][1], p = s[["Pr(>F)"]][1])
  }) |> list_rbind() |> save_tab("S5_type_comparison")

TYPE_COL <- set_names(c(PAL[["amber"]], PAL[["sage"]], PAL[["green"]]), CM_LAB)
TYPE_LAB <- set_names(c("Low citric/malic ratio", "Mid citric/malic ratio",
                        "High citric/malic ratio"), CM_LAB)
if (OPT[["ggtern"]]) {
  p_tern <- ggtern::ggtern(comp, ggtern::aes(x = malic, y = citric, z = oxalic)) +
    geom_point(aes(colour = 조성유형, size = total_acid), alpha = .8) +
    scale_colour_manual(values = TYPE_COL, name = "Citric/malic ratio tertile",
                        labels = TYPE_LAB) +
    labs(x = "Malic acid (%)", y = "Citric acid (%)", z = "Oxalic acid (%)")
  ggsave(file.path(PATH$fig, "S5_ternary.png"), p_tern, width = 8, height = 7,
         dpi = 300, bg = "white")
} else {
  # 구성비(%)끼리의 산점도는 합이 100으로 묶여 있어(compositional closure)
  # 음의 직선이 자동으로 나온다. 실제 상충 관계는 절대 함량에서 봐야 한다.
  message("[04] ggtern 미설치 -> 절대 함량 산점도로 대체")
  p_tern <- type_dat |>
    ggplot(aes(malic_mean, citric_mean)) +
    geom_point(aes(colour = 조성유형, size = oxalic_mean), alpha = .82) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey40", linewidth = .5, linetype = "dashed") +
    scale_colour_manual(values = TYPE_COL, name = "Citric/malic ratio tertile",
                        labels = TYPE_LAB) +
    scale_size_continuous(range = c(1, 5.5), name = OA_LAB_OXALIC) +
    labs(x = OA_LAB_MALIC, y = OA_LAB_CITRIC)
  save_fig(p_tern, "S5_composition_scatter", w = 8.5, h = 6.5)
}

# 조성비는 합이 100으로 고정된 폐쇄 자료다. 구성비끼리의 상관은 이 제약 때문에
# 자동으로 음이 되므로, 상충 관계의 근거로는 절대 함량 상관(S-7)을 쓴다.
comp |>
  summarise(across(c(oxalic, malic, citric), list(mean = mean, sd = sd)),
            r_pct_malic_citric  = cor(malic, citric),
            .by = 조성유형) |>
  mutate(r_abs_malic_citric = cor(type_dat$malic_mean, type_dat$citric_mean)) |>
  save_tab("S5_composition_summary")

## ---- S-6 표현형 그룹 간 비교 ------------------------------------------------
# 자원 하나 = 반복 하나. S-3(개체 반복)과는 실험 단위가 다른 별개 검정.
GROUP_VARS  <- c("fruit_shape_mode", "flower_diameter_mode", "leaf_green_intensity_mode",
                 "flesh_color_mode", "ground_color_mode")
TEST_VARS   <- c(ACID_MEAN, "total_acid")
MIN_GROUP_N <- 6

grp_dat <- oa_use |> mutate(across(all_of(GROUP_VARS), \(x) str_remove(x, "\\*$")))

grp_slice <- function(g, v) {
  grp_dat |>
    select(g = all_of(g), y = all_of(v)) |>
    filter(!is.na(g), is.finite(y)) |>
    filter(n() >= MIN_GROUP_N, .by = g) |>
    mutate(g = factor(g))
}

group_test <- expand_grid(그룹변수 = GROUP_VARS, 변수 = TEST_VARS) |>
  mutate(res = map2(그룹변수, 변수, \(g, v) {
    d <- grp_slice(g, v)
    if (n_distinct(d$g) < 2) return(NULL)
    s <- summary(aov(y ~ g, data = d))[[1]]
    tibble(k = n_distinct(d$g), n = nrow(d), F = s[["F value"]][1],
           p_anova = s[["Pr(>F)"]][1], KW_p = kruskal.test(y ~ g, data = d)$p.value,
           eta2 = s[["Sum Sq"]][1] / sum(s[["Sum Sq"]]))
  })) |>
  filter(!map_lgl(res, is.null)) |>
  unnest(res) |>
  mutate(p_adj = p.adjust(p_anova, "BH"), 유의 = p_adj < ALPHA) |>
  save_tab("S6_group_comparison")

message("[04] S-6 ANOVA: 유의 ", sum(group_test$유의), " / ", nrow(group_test), " 조합")
if (any(group_test$유의)) print(filter(group_test, 유의))

## 등가성 검정 (TOST) — 비유의를 "차이 없음"으로 쓰기 위한 절차
# 등가 한계 = Cohen's d 0.5 에 해당하는 원단위 폭
tost_pair <- function(x, y, eps_d = 0.5) {
  s_p <- sqrt(((length(x)-1)*var(x) + (length(y)-1)*var(y)) / (length(x)+length(y)-2))
  eps <- eps_d * s_p
  tibble(diff = mean(x) - mean(y), eps = eps,
         tost_p = max(t.test(x, y, mu = -eps, alternative = "greater")$p.value,
                      t.test(x, y, mu =  eps, alternative = "less")$p.value))
}

equiv_tab <- expand_grid(그룹변수 = GROUP_VARS, 변수 = TEST_VARS) |>
  mutate(res = map2(그룹변수, 변수, \(g, v) {
    d  <- grp_slice(g, v)
    lv <- levels(droplevels(d$g))
    if (length(lv) < 2) return(NULL)
    combn(lv, 2, simplify = FALSE) |>
      map(\(pr) tost_pair(d$y[d$g == pr[1]], d$y[d$g == pr[2]]) |>
            mutate(비교 = str_c(pr, collapse = " vs "))) |>
      list_rbind()
  })) |>
  filter(!map_lgl(res, is.null)) |>
  unnest(res) |>
  mutate(tost_p_adj = p.adjust(tost_p, "BH"),
         등가 = tost_p < ALPHA, 등가_BH = tost_p_adj < ALPHA) |>
  save_tab("S6_equivalence_TOST")

message("[04] S-6 TOST: 등가 판정 ", sum(equiv_tab$등가), " / ", nrow(equiv_tab), " 쌍",
        " (BH 보정 후 ", sum(equiv_tab$등가_BH), "쌍)")

# ANOVA 와 TOST 를 결합한 쌍별 결론 — 콘솔 메시지로만 두지 않고 표로 남긴다.
conclusion_tab <- equiv_tab |>
  left_join(select(group_test, 그룹변수, 변수, p_anova, p_adj_anova = p_adj),
            by = c("그룹변수", "변수")) |>
  mutate(결론 = case_when(
    p_adj_anova < ALPHA ~ "차이 있음",
    등가                ~ "실질적으로 차이 없음 (등가)",
    .default            = "판정 불가 (표본 부족)")) |>
  select(그룹변수, 변수, 비교, diff, eps, p_anova, p_adj_anova, tost_p, tost_p_adj, 결론) |>
  save_tab("S6_pairwise_conclusion")

message("     쌍별 결론: ", str_c(names(table(conclusion_tab$결론)), " ",
                                  as.integer(table(conclusion_tab$결론)), collapse = " | "))
message("     ANOVA 비유의 + TOST 유의 = 실질적으로 차이 없음.")
message("     둘 다 비유의면 표본이 부족한 것이지 차이가 없는 것이 아니다.")

p_box <- grp_dat |>
  select(no, all_of(GROUP_VARS), all_of(ACID_MEAN)) |>
  pivot_longer(all_of(GROUP_VARS), names_to = "그룹변수", values_to = "수준") |>
  pivot_longer(all_of(ACID_MEAN), names_to = "acid", values_to = "value") |>
  filter(!is.na(수준)) |>
  mutate(acid = relabel_acid(acid)) |>
  ggplot(aes(수준, value)) +
  geom_boxplot(outlier.size = .7, fill = PAL[["pale"]], colour = PAL[["green"]], linewidth = .35) +
  facet_grid(acid ~ 그룹변수, scales = "free", labeller = labeller(acid = label_parsed)) +
  labs(x = NULL, y = OA_LAB) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
save_fig(p_box, "S6_group_boxplot", w = 13, h = 8)

## ---- S-7 형질 간 상관 -------------------------------------------------------
cor_vars <- c(ACID_MEAN, DERIVED, "ta_mean", "ssc_mean", "fruit_weight_mean",
              "bloom_doy_mean") |> keep(\(x) x %in% names(oa_use))

# BH 보정은 고유 쌍에만 적용한다. 대각선(자기상관 p = 0)과 대칭 중복까지 넣으면
# 검정 수가 부풀려져 보정 p 가 체계적으로 작아진다(반보수적).
cor_pairs <- t(combn(cor_vars, 2)) |> as_tibble(.name_repair = ~ c("x", "y")) |>
  mutate(res = map2(x, y, \(a, b) {
    ok <- is.finite(oa_use[[a]]) & is.finite(oa_use[[b]])
    if (sum(ok) < 10) return(tibble(n = sum(ok), pearson = NA_real_,
                                    spearman = NA_real_, p = NA_real_))
    ct <- suppressWarnings(cor.test(oa_use[[a]][ok], oa_use[[b]][ok]))
    tibble(n = sum(ok), pearson = unname(ct$estimate),
           spearman = suppressWarnings(cor(oa_use[[a]][ok], oa_use[[b]][ok],
                                           method = "spearman")),
           p = ct$p.value)
  })) |>
  unnest(res) |>
  mutate(p_adj = p.adjust(p, "BH"))          # 고유 36쌍 기준

message("[04] S-7 BH 보정 대상 고유 쌍 ", nrow(cor_pairs), "개")

# 히트맵·표는 대칭 격자로 펼치되 보정 p 는 위에서 계산한 값을 그대로 복사한다.
cor_long <- bind_rows(
  cor_pairs,
  rename(cor_pairs, x = y, y = x),
  tibble(x = cor_vars, y = cor_vars, n = map_int(cor_vars, \(v) sum(is.finite(oa_use[[v]]))),
         pearson = 1, spearman = 1, p = NA_real_, p_adj = NA_real_)
) |>
  arrange(match(x, cor_vars), match(y, cor_vars)) |>
  save_tab("S7_correlation")

p_cor <- ggplot(cor_long, aes(x, y, fill = pearson)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = if_else(x == y, "",
                                sprintf("%.2f%s", pearson, if_else(p_adj < .05, "*", "")))),
            size = 2.8, colour = "grey12") +
  scale_fill_gradient2(low = PAL[["rust"]], mid = "grey96", high = PAL[["green"]],
                       midpoint = 0, limits = c(-1, 1),
                       name = "Pearson r\n(* BH-adjusted p < 0.05)") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
save_fig(p_cor, "S7_correlation_heatmap", w = 8.5, h = 7.5)

ta_check <- cor_long |> filter(x == "total_acid", y == "ta_mean")
message(sprintf("[04] S-7 총산 vs titratable_acidity: r = %.3f (p = %.3f, n = %d)",
                ta_check$pearson, ta_check$p, ta_check$n))
message("     양(+) 상관이 아니면 두 측정의 시료·단위 정합성을 확인할 것.")

## ---- S-8 목적별 선발 --------------------------------------------------------
sel_base <- ci_tab |>
  select(acid, no, mean, lwr, upr) |>
  pivot_wider(names_from = acid, values_from = c(mean, lwr, upr),
              names_glue = "{acid}_{.value}") |>
  left_join(select(oa_use, no, total_acid, citric_malic_ratio, qc_flag), by = "no") |>
  left_join(select(comp, no, 조성유형), by = "no")

# 총산은 세 산의 합이므로 자체 오차항이 없다. 세 산의 풀링 MSE 를 더해 표준오차를
# 전파한다: Var(총산 평균) = sum(MSE_i) / n_rep (세 측정이 독립이라는 가정).
# 이 열이 없으면 아래 pick() 의 이름 치환이 실패해 CI 가 [값, 값] 으로 퇴화한다.
se_total <- sqrt(sum(anova_tab$MS_within) / N_REP)
t_total  <- qt(1 - ALPHA/2, min(anova_tab$df_within))
sel_base <- sel_base |>
  mutate(total_acid_lwr = total_acid - t_total * se_total,
         total_acid_upr = total_acid + t_total * se_total)

pick <- function(col, n = 10, desc = TRUE, label) {
  lwr_col <- str_c(str_remove(col, "_mean$"), "_lwr")
  upr_col <- str_c(str_remove(col, "_mean$"), "_upr")
  stopifnot(all(c(lwr_col, upr_col) %in% names(sel_base)))
  # mg/g 단위에서 옥살산은 0.2~1.0 범위다. 소수 1자리로는 등급이 뭉개지므로
  # 자릿수를 값의 크기에 맞춘다.
  dg <- if (max(sel_base[[col]], na.rm = TRUE) < 5) 3L else 2L
  d <- arrange(sel_base, no)
  d <- if (desc) slice_max(d, .data[[col]], n = n, with_ties = FALSE)
       else      slice_min(d, .data[[col]], n = n, with_ties = FALSE)
  d |> transmute(목표 = label, no, 값 = round(.data[[col]], dg),
                 CI = str_glue("[{round(.data[[lwr_col]], dg)}, {round(.data[[upr_col]], dg)}]"),
                 조성유형, qc_flag,
                 재검증권장 = str_detect(replace_na(qc_flag, ""), "high_rsd"))
}

selection_oa <- list(
  pick("oxalic_acid_mean", 10, FALSE, "저옥살산 (항영양인자 저감)"),
  pick("citric_acid_mean", 10, TRUE,  "고구연산 (청량 산미)"),
  pick("total_acid",       10, TRUE,  "고총산 (매실청·매실주 원료)"),
  pick("total_acid",       10, FALSE, "저총산 (생식·저산 가공)")
) |> list_rbind() |> save_tab("S8_selection_organic_acid")

message("[04] S-8 선발 ", nrow(selection_oa), "행 | 재검증 권장(high_rsd) ",
        sum(selection_oa$재검증권장), "건")
message("[04] S-1 ~ S-8 완료")
