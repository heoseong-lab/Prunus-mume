# =============================================================================
# 01_preprocess.R — 계획서 P-2 ~ P-6
#   P-2 범주 정제·순서형 코딩 / P-3 날짜 파생 / P-4 이상값·결측
#   P-5 연차 효과 보정(BLUP) / P-6 자원 단위 분석 매트릭스
# 산출: dat_long (484행 정제) · acc_mat (204 x 형질) · repeatability_tab
# =============================================================================

if (!exists("PATH")) source(file.path("R", "00_setup.R"), encoding = "UTF-8")

raw <- read_csv(PATH$raw, show_col_types = FALSE,
                locale = locale(encoding = "UTF-8"), na = c("", "NA")) |>
  mutate(across(everything(), \(x) na_if(str_trim(x), "")))

message("[01] 원자료 ", nrow(raw), "행 x ", ncol(raw), "열")

## ---- P-2 범주 정제 ----------------------------------------------------------
# 병해 등급: 연차마다 구간이 다르다 (2026년만 "1~10%", 그 외 "1~5%"/"5~10%";
# "10~20%"와 "10~50%" 중복). 공통 4단계로 통합.
harmonize_severity <- function(x) {
  case_when(
    is.na(x)                             ~ NA_character_,
    x == "0%"                            ~ "0%",
    x == "~1%"                           ~ "~1%",
    x %in% c("1~5%", "5~10%", "1~10%")   ~ "1~10%",
    x %in% c("10~20%", "10~50%", "20%~") ~ "10%~",
    .default = NA_character_
  )
}

## ---- P-2a 단위 보정 ---------------------------------------------------------
# 2024년 fruit_weight 는 과실 1개가 아니라 3개를 합쳐 단 무게로 기록되어 있다
# (조사자 확인, 2026-08-27). 개당 무게로 환산해야 연차 간 비교가 성립한다.
FW_DIVISOR_2024 <- 3

dat_long <- raw |>
  mutate(
    year = as.integer(year),
    fruit_weight = if_else(year == 2024L, as.numeric(fruit_weight) / FW_DIVISOR_2024,
                           as.numeric(fruit_weight)),
    skin_color_intensity = case_when(                               # 철자 오류
      skin_color_intensity == "stong"      ~ "strong",
      skin_color_intensity == "very stong" ~ "very strong",
      .default = skin_color_intensity),
    leaf_green_intensity = na_if(leaf_green_intensity, "other"),    # 순서형에 못 넣음
    brown_rot_severity   = harmonize_severity(brown_rot_severity),
    shot_hole_severity   = harmonize_severity(shot_hole_severity),
    across(all_of(trait_of("continuous")), as.numeric)
  )

# 극소 빈도 수준 병합 (다양성 지수 왜곡 방지)
# 규칙
#   순서형: 관측 수가 가장 적은 희소 수준부터 '가장 가까운 비희소 수준'으로 병합하고,
#           희소 수준이 남지 않을 때까지 반복한다. 인접 수준이 함께 희소이면
#           한 번에 병합되지 않으므로 반드시 반복이 필요하다.
#   명목형: 순서가 없어 인접 병합이 불가능하므로 희소 수준을 "other" 로 통합하되,
#           통합해도 MIN_LEVEL_N 에 못 미치면(희소 수준이 하나뿐인 경우) 통합이
#           수준 수를 줄이지 못하고 이름만 가리므로 원 수준을 그대로 둔다.
MIN_LEVEL_N <- 5
merge_log <- list()
for (tr in c(trait_of("nominal"), trait_of("ordinal"))) {
  is_ord <- TRAITS$type[match(tr, TRAITS$trait)] == "ordinal"
  if (is_ord) {
    lv <- levels_of(tr)
    repeat {
      tb   <- table(factor(dat_long[[tr]], levels = lv))
      rare <- names(tb)[tb > 0 & tb < MIN_LEVEL_N]
      ok   <- names(tb)[tb >= MIN_LEVEL_N]
      if (!length(rare) || !length(ok)) break
      r   <- rare[which.min(tb[rare])]                       # 가장 드문 수준부터
      tgt <- ok[which.min(abs(match(ok, lv) - match(r, lv)))] # 가장 가까운 비희소 수준
      dat_long[[tr]][dat_long[[tr]] == r] <- tgt
      merge_log[[tr]] <- c(merge_log[[tr]], str_c(r, " -> ", tgt))
    }
  } else {
    tb   <- table(dat_long[[tr]])
    rare <- names(tb)[tb < MIN_LEVEL_N]
    if (length(rare) && sum(tb[rare]) >= MIN_LEVEL_N) {
      dat_long[[tr]][dat_long[[tr]] %in% rare] <- "other"
      merge_log[[tr]] <- c(merge_log[[tr]], str_c(str_c(rare, collapse = "+"), " -> other"))
    } else if (length(rare)) {
      merge_log[[tr]] <- c(merge_log[[tr]],
                           str_c(str_c(rare, collapse = "+"), " -> 유지 (통합해도 n < ",
                                 MIN_LEVEL_N, ")"))
    }
  }
}
if (length(merge_log)) {
  message("[01] P-2 희소수준 병합:")
  iwalk(merge_log, \(v, k) message("     ", k, ": ", str_c(v, collapse = " | ")))
}
imap(merge_log, \(v, k) tibble(trait = k, 처리 = v)) |>
  list_rbind() |>
  save_tab("P2_level_merge_log")

# 병합 후 잔존 희소 수준 점검 (감사용)
walk(c(trait_of("nominal"), trait_of("ordinal")), \(tr) {
  tb <- table(dat_long[[tr]])
  if (any(tb < MIN_LEVEL_N))
    message("     [잔존] ", tr, ": ",
            str_c(names(tb)[tb < MIN_LEVEL_N], "=", tb[tb < MIN_LEVEL_N], collapse = ", "))
})

## ---- P-3 날짜 -> 수치 -------------------------------------------------------
to_doy <- function(x, max_month = 5L) {
  d <- suppressWarnings(as.Date(x))
  d[!is.na(d) & month(d) > max_month] <- NA      # KJM101 2025-06-25 등 이상값
  yday(d)
}

dat_long <- dat_long |>
  mutate(
    bloom_doy     = to_doy(full_bloom_date, 5L),
    maturity_doy  = to_doy(maturity_date,   8L),
    dev_days      = maturity_doy - bloom_doy,
    # 성숙일은 자원별 성숙기가 아니라 연 4~5회 일괄 수확일 -> 배치 공변량으로만
    harvest_batch = factor(str_c(year, "_", maturity_date)),
    accession     = factor(no),
    year_f        = factor(year)
  )

message("[01] 만개일 이상값 제거 ",
        sum(!is.na(dat_long$full_bloom_date) & is.na(dat_long$bloom_doy)),
        "행 | 수확 배치 ", nlevels(dat_long$harvest_batch), "개")

## ---- P-4 이상값 스크리닝 ----------------------------------------------------
# 자동 제거하지 않는다. 원장 대조용 목록만 만든다.
c(trait_of("continuous"), "bloom_doy", "dev_days") |>
  map(\(tr) dat_long |>
        filter(flag_outlier(.data[[tr]])) |>
        transmute(trait = tr, no, year, value = .data[[tr]])) |>
  list_rbind() |>
  save_tab("P4_outlier_report")

# fruit_weight 는 보정 전 2024년 값이 3배로 기록되어 있었고(33.2 g), 보정 후에도
# 연차 평균이 11.05 / 10.53 / 9.88 g 로 조금씩 이동한다. 전체 분포 기준 스크리닝은
# 이 이동에 끌려가므로 연차 내 기준으로 다시 본다.
dat_long |>
  group_by(year) |>
  filter(flag_outlier(fruit_weight)) |>
  ungroup() |>
  select(no, year, fruit_weight) |>
  save_tab("P4_fruit_weight_outlier_within_year")

## ---- P-5 연차 효과 보정 (BLUP) ----------------------------------------------
# y = mu + Y(연차, fixed) + G(자원, random) + B(수확배치, random) + e
fit_blup <- function(trait, data, use_batch = TRUE) {
  d <- data |>
    select(accession, year_f, harvest_batch, y = all_of(trait)) |>
    filter(is.finite(y))
  if (n_distinct(d$accession) < 10) return(NULL)

  form <- if (use_batch && n_distinct(d$harvest_batch) > 2) {
    y ~ year_f + (1 | accession) + (1 | harvest_batch)
  } else {
    y ~ year_f + (1 | accession)
  }
  fit <- try(lmer(form, data = d,
                  control = lmerControl(check.conv.singular = "ignore")), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  vc    <- as_tibble(VarCorr(fit))
  var_g <- vc |> filter(grp == "accession") |> pull(vcov)
  var_e <- vc |> filter(grp == "Residual")  |> pull(vcov)
  var_b <- vc |> filter(grp == "harvest_batch") |> pull(vcov)
  if (!length(var_b)) var_b <- 0

  re   <- ranef(fit, condVar = TRUE)$accession
  blup <- re |>
    rownames_to_column("no") |>
    as_tibble() |>
    transmute(no, "{trait}_blup" := `(Intercept)` + fixef(fit)[[1]])
  # 자원별 BLUP 신뢰도 1 - PEV_i/σ²G. Spearman-Brown(n̄) 이 실제 예측값의 신뢰도를
  # 얼마나 잘 근사하는지 대조하기 위해 평균을 함께 보고한다 (Cullis et al. 2006 의 정신).
  pev <- attr(re, "postVar")[1, 1, ]
  rel_pev <- mean(1 - pev / var_g)

  # 반복성의 정의 — 수확 배치는 자원의 성질이 아니라 조사 일정의 인공물이므로
  # 주 지표는 배치 분산을 제거한 '조정 반복성' R = σ²G / (σ²G + σ²e) 로 둔다.
  # 다만 같은 자원의 연차 간 측정은 서로 다른 배치에서 이루어지므로 배치 분산도
  # 비유전 변동의 일부다. 정의에 따라 값이 크게 달라지는 형질이 있어(과중 0.369 vs
  # 0.219) 배치를 분모에 포함한 보수적 정의도 함께 보고한다.
  lst(trait, blup, var_g, var_b, var_e, rel_pev,
      repeatability = var_g / (var_g + var_e),
      repeatability_incl_batch = var_g / (var_g + var_b + var_e),
      n_obs = nrow(d), n_acc = n_distinct(d$accession))
}

# 개화기·과실발육일수도 BLUP 대상에 포함한다. 단 dev_days 는 성숙일에서 파생되므로
# 수확배치를 임의효과로 넣으면 순환이 되어 제외한다.
BLUP_TRAITS <- c(trait_of("continuous"), "bloom_doy", "dev_days")

blup_fits <- BLUP_TRAITS |>
  map(\(tr) fit_blup(tr, dat_long, use_batch = tr %in% trait_of("continuous"))) |>
  compact()
blup_fits <- set_names(blup_fits, map_chr(blup_fits, "trait"))

repeatability_tab <- blup_fits |>
  map(\(f) tibble(trait = f$trait, n_obs = f$n_obs, n_accession = f$n_acc,
                  n_bar = f$n_obs / f$n_acc,
                  var_G = f$var_g, var_batch = f$var_b, var_e = f$var_e,
                  repeatability = f$repeatability,
                  repeatability_incl_batch = f$repeatability_incl_batch,
                  blup_reliability_pev = f$rel_pev)) |>
  list_rbind() |>
  arrange(desc(repeatability)) |>
  mutate(판정 = case_when(repeatability >= .5 ~ "채택",
                          repeatability >= .2 ~ "가중치 축소",
                          .default = "제외 검토")) |>
  save_tab("A2_repeatability")

print(repeatability_tab)

## ---- P-6 자원 단위 매트릭스 -------------------------------------------------
# 명목형: 최빈값(+동점 플래그) / 순서형: UPOV 코드 중앙값 / 연속형: BLUP
cat_summary <- dat_long |>
  summarise(across(all_of(trait_of("nominal")),
                   list(mode = mode_chr, tie = is_tied), .names = "{.col}__{.fn}"),
            across(all_of(trait_of("ordinal")), mode_chr),
            .by = no)

ord_summary <- dat_long |>
  select(no, all_of(trait_of("ordinal"))) |>
  pivot_longer(-no, names_to = "trait", values_to = "level") |>
  filter(!is.na(level)) |>
  mutate(code = map2_dbl(level, trait, \(l, t) upov_code(l, levels_of(t)))) |>
  summarise(code = median(code, na.rm = TRUE), .by = c(no, trait)) |>
  pivot_wider(names_from = trait, values_from = code, names_glue = "{trait}_code")

acc_base <- dat_long |>
  summarise(n_years   = n_distinct(year),
            years     = str_c(sort(unique(year)), collapse = ";"),
            # 연차 미보정 단순 평균 — BLUP(_blup) 과 구분되도록 이름에 표시한다.
            bloom_doy_rawmean = mean(bloom_doy, na.rm = TRUE),
            dev_days_rawmean  = mean(dev_days,  na.rm = TRUE),
            .by = no) |>
  mutate(across(c(bloom_doy_rawmean, dev_days_rawmean),
                \(x) if_else(is.nan(x), NA_real_, x)),
         single_year = n_years == 1) |>
  left_join(cat_summary, by = "no") |>
  left_join(ord_summary, by = "no")

acc_mat <- reduce(map(blup_fits, "blup"), left_join, by = "no", .init = acc_base)

save_tab(acc_mat, "P6_analysis_matrix")
message("[01] 자원 매트릭스 ", nrow(acc_mat), "행 x ", ncol(acc_mat), "열 | 단년 자료 ",
        sum(acc_mat$single_year), "점")
