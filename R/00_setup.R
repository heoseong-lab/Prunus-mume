# =============================================================================
# 00_setup.R — 매실 유전자원 분석: 공통 설정 · 패키지 · 헬퍼 함수
# =============================================================================
# 실행: source("R/run_all.R", encoding = "UTF-8")  또는  Rscript R/run_all.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)      # dplyr tidyr purrr readr stringr forcats tibble ggplot2
  library(lme4)           # 혼합모형 · BLUP
  library(car)            # 등분산 검정
  library(agricolae)      # LSD · Duncan · Tukey 평균분리
  library(FactoMineR)     # FAMD
  library(factoextra)     # 다변량 시각화
  library(cluster)        # Gower 거리 · silhouette
  library(ggrepel)
})
# scales 는 attach 하지 않는다 — scales::discard 가 purrr::discard 를 가린다.
# 필요한 곳에서 scales:: 로 호출.

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# 선택 패키지 — 없으면 대체 경로로 진행
OPT <- c(clustMixType = has_pkg("clustMixType"),   # k-prototypes 교차검증 (A-5)
         ggtern       = has_pkg("ggtern"),         # 삼각도 (S-5)
         ordinal      = has_pkg("ordinal"))        # 순서형 CLMM 연차 효과 검정 (A-2)
if (any(!OPT)) {
  message("[00] 선택 패키지 미설치: ", paste(names(OPT)[!OPT], collapse = ", "),
          "  -> 대체 경로로 진행. 설치하려면: install.packages(c(",
          paste0('"', names(OPT)[!OPT], '"', collapse = ", "), "))")
}

## ---- 경로 -------------------------------------------------------------------
find_root <- function() {
  for (cand in c(getwd(), dirname(getwd()), file.path(getwd(), ".."))) {
    if (file.exists(file.path(cand, "매실 전체 데이터.csv"))) return(normalizePath(cand))
  }
  stop("프로젝트 루트를 찾지 못했습니다. 작업 디렉터리를 '매실 유전자원' 폴더로 설정하세요.")
}
PROJ_ROOT <- find_root()

PATH <- list(
  raw    = file.path(PROJ_ROOT, "매실 전체 데이터.csv"),
  oa     = file.path(PROJ_ROOT, "매실_유기산_자원.csv"),
  oa_raw = file.path(PROJ_ROOT, "유기산_반복원자료.csv"),   # 있으면 우선 사용
  fig    = file.path(PROJ_ROOT, "output", "figures"),
  tab    = file.path(PROJ_ROOT, "output", "tables")
)
walk(PATH[c("fig", "tab")], dir.create, recursive = TRUE, showWarnings = FALSE)

## ---- 테마 -------------------------------------------------------------------
PAL <- c(green = "#2F5D3A", amber = "#7E560F", rust = "#94351F",
         sage = "#5C6B57", pale = "#DBE7D4")

# 그림에는 제목·부제를 넣지 않는다. 저널 조판에서 그 역할은 본문 캡션이 맡으며,
# 그림 안에 제목이 들어가면 캡션과 중복된다.
theme_set(
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_blank(),
          plot.subtitle    = element_blank(),
          strip.text       = element_text(face = "bold"))
)
set.seed(2026)

## ---- 형질 사전 (P-1) --------------------------------------------------------
# type: continuous / ordinal / nominal / date / id
# levels: 순서형은 낮은 등급 -> 높은 등급 순 (UPOV 1·3·5·7·9 코딩에 사용)
TRAITS <- tribble(
  ~trait,                    ~type,         ~levels,
  "no",                      "id",          NULL,
  "year",                    "id",          NULL,
  "full_bloom_date",         "date",        NULL,
  "maturity_date",           "date",        NULL,
  "fruit_shape",             "nominal",     c("circular","oblate","elliptic","ovate","obovate","oblong"),
  "ground_color",            "nominal",     c("light green","medium green","yellow","red"),
  "skin_color_intensity",    "ordinal",     c("absent or very weak","weak","medium","strong","very strong"),
  "flesh_color",             "nominal",     c("whitish green","white","cream","light orange","orange","dark orange"),
  "stone_adherence",         "ordinal",     c("clingstone","semi-clingstone","semi-freestone","freestone"),
  "leaf_green_intensity",    "ordinal",     c("light","medium","dark"),
  "leaf_base_shape",         "nominal",     c("acute","obtuse","truncate","cordate"),
  "leaf_pubescence",         "ordinal",     c("absent","weak","medium","strong"),
  "leaf_angle",              "ordinal",     c("acute","right-angled","moderately obtuse","strongly obtuse"),
  "leaf_length",             "ordinal",     c("short","medium","long"),
  "leaf_width",              "ordinal",     c("narrow","medium","broad"),
  "fruit_weight",            "continuous",  NULL,
  "flower_diameter",         "ordinal",     c("small","medium","large"),
  "petal_number",            "continuous",  NULL,
  # 병해 2형질은 연차별 등급 체계가 달라 공통 4단계로 재코딩(01_preprocess)
  "brown_rot_severity",      "ordinal",     c("0%","~1%","1~10%","10%~"),
  "shot_hole_severity",      "ordinal",     c("0%","~1%","1~10%","10%~"),
  "soluble_solids_content",  "continuous",  NULL,
  "titratable_acidity",      "continuous",  NULL
)

trait_of  <- function(ty) TRAITS |> filter(type == ty) |> pull(trait)
levels_of <- function(tr) TRAITS$levels[[match(tr, TRAITS$trait)]]

## ---- 헬퍼 -------------------------------------------------------------------
# 최빈값 (동점이면 알파벳 순 첫 값)
mode_chr <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}
is_tied <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) < 2) return(FALSE)
  tb <- sort(table(x), decreasing = TRUE)
  length(tb) > 1 && tb[[1]] == tb[[2]]
}

# 순서형 -> UPOV 1·3·5·7·9 (수준 수에 맞춰 균등 배치)
upov_code <- function(x, lv) {
  i <- match(as.character(x), lv)
  k <- length(lv)
  if (k <= 1) return(rep(5, length(i)))
  round(1 + (i - 1) * 8 / (k - 1))
}

# Shannon–Weaver 다양성지수 (0~1 정규화). 연속형은 구간화 후 계산.
# 정규화 분모 k 는 '정의된 등급 수'다 — 연속형은 구간 수(bins), 범주형은 병합 후 등급 수.
# (관측된 비영 계급 수로 나누면 빈 구간이 있는 형질의 H' 가 과대평가된다)
shannon <- function(x, bins = 10) {
  if (is.numeric(x)) {
    x <- x[is.finite(x)]
    if (n_distinct(x) < 2) return(0)
    x <- cut(x, breaks = bins)
    k <- bins
  } else {
    k <- n_distinct(x[!is.na(x)])
  }
  p <- table(x[!is.na(x)]) |> prop.table()
  p <- p[p > 0]
  if (length(p) < 2 || k < 2) return(0)
  as.numeric(-sum(p * log(p)) / log(k))
}

# Fleiss' kappa — 연차를 평가자로 본 다중평가자 일치도 (우연 보정)
# df: 한 형질에 대해 subject(자원) x rater(연차) 형태의 long 자료 (id, level)
# 관측 수가 동일한 자원만 사용한다.
fleiss_kappa <- function(id, level, min_raters = 2L) {
  d  <- tibble(id, level) |> filter(!is.na(level))
  m  <- d |> count(id) |> filter(n == max(n)) |> pull(id)   # 최다 관측 수 자원만
  d  <- d |> filter(id %in% m)
  R  <- d |> count(id) |> pull(n) |> first()
  if (is.na(R) || R < min_raters || n_distinct(d$level) < 2) return(NA_real_)

  tab <- table(d$id, d$level)                  # n_ij
  N   <- nrow(tab)
  P_i <- (rowSums(tab^2) - R) / (R * (R - 1))
  P_bar <- mean(P_i)
  P_e   <- sum((colSums(tab) / (N * R))^2)
  if (P_e >= 1) return(NA_real_)
  as.numeric((P_bar - P_e) / (1 - P_e))
}

# κ 산출에 실제로 사용된 자원 수 — 위 함수가 최다 관측 자원만 쓰므로 별도로 보고한다.
fleiss_kappa_n <- function(id, level) {
  d <- tibble(id, level) |> filter(!is.na(level))
  if (!nrow(d)) return(NA_integer_)
  cnt <- count(d, id)
  as.integer(sum(cnt$n == max(cnt$n)))
}

# Krippendorff's α — 자료형에 맞는 불일치 함수(nominal / ordinal / interval)를 쓰는 우연 보정 일치도.
# Fleiss κ 와의 차이: (1) 순서형의 순서를 반영한다(1등급 차이와 4등급 차이를 구분).
# (2) 자원마다 관측 수가 달라도 2회 이상 관측된 자원을 모두 쓴다(κ 는 최다 관측 자원만 씀).
# (3) 연속형에도 같은 틀(interval)을 적용할 수 있다. Krippendorff (2019) Content Analysis 4판 12장.
# value: nominal 은 임의 정수 코드, ordinal 은 낮은 등급 -> 높은 등급 순의 정수 코드, interval 은 수치.
kripp_alpha <- function(id, value, level = c("nominal", "ordinal", "interval")) {
  level <- match.arg(level)
  d <- tibble(id, value) |> filter(!is.na(value)) |> add_count(id, name = "m") |> filter(m >= 2)
  if (!nrow(d) || n_distinct(d$value) < 2) return(NA_real_)
  vals <- sort(unique(d$value)); K <- length(vals)
  o <- matrix(0, K, K)                                   # 관측 일치행렬 (unit 별 쌍, 1/(m-1) 가중)
  for (u in split(d$value, d$id)) {
    ix <- match(u, vals); mu <- length(u)
    for (a in seq_len(mu)) for (b in seq_len(mu)) if (a != b)
      o[ix[a], ix[b]] <- o[ix[a], ix[b]] + 1 / (mu - 1)
  }
  nc <- rowSums(o); n <- sum(nc)
  d2 <- matrix(0, K, K)                                  # 불일치 함수 δ²
  for (i in seq_len(K)) for (j in seq_len(K)) d2[i, j] <- switch(level,
    nominal  = as.numeric(i != j),
    interval = (vals[i] - vals[j])^2,
    ordinal  = if (i == j) 0 else { lo <- min(i, j); hi <- max(i, j)
                                    (sum(nc[lo:hi]) - (nc[lo] + nc[hi]) / 2)^2 })
  Do <- sum(o * d2); De <- sum(outer(nc, nc) * d2) / (n - 1)
  if (!is.finite(De) || De <= 0) return(NA_real_)
  as.numeric(1 - Do / De)
}
kripp_alpha_n <- function(id, value) {
  d <- tibble(id, value) |> filter(!is.na(value)) |> count(id) |> filter(n >= 2)
  nrow(d)
}

# Spearman-Brown: 단일 관측 신뢰도 R 을 n 회 평균값의 신뢰도로 환산.
# 일원 확률효과 모형에서 자원별 BLUP 신뢰도 1 - PEV_i/σ²G 는 정확히 이 식과 같다(n = n_i).
# 음수 입력은 그대로 돌려준다 — 가중치는 어차피 0 이 되고, 표에는 원래 부호를 남기기 위함.
spearman_brown <- function(R, n) if_else(R > 0, n * R / (1 + (n - 1) * R), R)

# Adjusted Rand Index — 두 군집 결과의 일치도 (우연 보정)
adj_rand <- function(a, b) {
  tab <- table(a, b)
  ch2 <- function(x) x * (x - 1) / 2
  idx <- sum(ch2(tab))
  ea  <- sum(ch2(rowSums(tab))); eb <- sum(ch2(colSums(tab)))
  exp <- ea * eb / ch2(sum(tab))
  (idx - exp) / (0.5 * (ea + eb) - exp)
}

cramers_v <- function(x, y) {
  tb <- table(x, y)
  if (any(dim(tb) < 2)) return(NA_real_)
  # V 는 무보정 χ² 기준이 표준이다 (2×2 에서 Yates 보정은 V 를 하향 편향시킨다)
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tb) * (min(dim(tb)) - 1))))
}

# correlation ratio eta^2 (연속 x 범주)
eta_squared <- function(num, grp) {
  keep <- is.finite(num) & !is.na(grp)
  num <- num[keep]; grp <- droplevels(factor(grp[keep]))
  if (nlevels(grp) < 2 || length(num) < 3) return(NA_real_)
  ss_t <- sum((num - mean(num))^2)
  if (ss_t == 0) return(NA_real_)
  ss_b <- tapply(num, grp, \(v) length(v) * (mean(v) - mean(num))^2) |> sum()
  as.numeric(ss_b / ss_t)
}

# Cramér's V 편향보정 (Bergsma 2013)
# V 는 범주 수가 크고 표본이 작을수록 위로 편향되어 독립일 때에도 기댓값이 0 보다 크다.
cramers_v_bc <- function(x, y) {
  tb <- table(x, y)
  if (any(dim(tb) < 2)) return(NA_real_)
  n <- sum(tb); r <- nrow(tb); k <- ncol(tb)
  if (n <= 1) return(NA_real_)
  chi  <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  phi2 <- as.numeric(chi) / n
  phi2_t <- max(0, phi2 - (r - 1) * (k - 1) / (n - 1))   # 귀무 기댓값만큼 뺀다
  r_t <- r - (r - 1)^2 / (n - 1)
  k_t <- k - (k - 1)^2 / (n - 1)
  den <- min(r_t, k_t) - 1
  if (den <= 0) return(NA_real_)
  sqrt(phi2_t / den)
}

# 조정 eta^2 — 집단 수 k 가 크면 eta^2 도 위로 편향된다. 조정 R^2 와 같은 형태로 보정.
eta_squared_adj <- function(num, grp) {
  keep <- is.finite(num) & !is.na(grp)
  num <- num[keep]; grp <- droplevels(factor(grp[keep]))
  n <- length(num); k <- nlevels(grp)
  if (k < 2 || n < 3 || n - k <= 0) return(NA_real_)
  e2 <- eta_squared(num, grp)
  if (is.na(e2)) return(NA_real_)
  max(0, 1 - (1 - e2) * (n - 1) / (n - k))
}

flag_outlier <- function(x, k = 3) {
  q <- quantile(x, c(.25, .75), na.rm = TRUE)
  iqr <- q[[2]] - q[[1]]
  !is.na(x) & (x < q[[1]] - k * iqr | x > q[[2]] + k * iqr)
}

# 모집단 표준편차 기준 표준화 잔차의 3·4차 적률
std_moment <- function(x, p) {
  x <- x[is.finite(x)]
  z <- (x - mean(x)) / (sd(x) * sqrt((length(x) - 1) / length(x)))
  mean(z^p)
}
skewness <- function(x) std_moment(x, 3)
kurtosis <- function(x) std_moment(x, 4) - 3

save_tab <- function(df, name) {
  write_excel_csv(df, file.path(PATH$tab, str_c(name, ".csv")))   # UTF-8 BOM (Excel 한글)
  invisible(df)
}
# 저장 직전에 제목·부제·캡션을 한 번 더 지운다. factoextra(fviz_*)·ggtern 은
# 호출하지 않아도 자체 기본 제목을 붙이므로 여기서 일괄 차단한다.
save_fig <- function(plot, name, w = 8, h = 6, dpi = 300) {
  if (inherits(plot, "ggplot"))
    plot <- plot + labs(title = NULL, subtitle = NULL, caption = NULL)
  ggsave(file.path(PATH$fig, str_c(name, ".png")), plot,
         width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(plot)
}

## ---- 그림용 형질 표기 (영문 저널 제출용) -------------------------------------
# 그림 축에 변수명(flower_diameter, fruit_weight_blup) 대신 원고 표기를 쓴다.
# "_blup" 접미사는 떼고 찾으며, 사전에 없는 이름은 그대로 돌려준다.
TRAIT_LABELS <- c(
  flower_diameter        = "Flower diameter",
  leaf_pubescence        = "Leaf pubescence",
  leaf_angle             = "Leaf angle",
  leaf_width             = "Leaf width",
  leaf_length            = "Leaf length",
  leaf_green_intensity   = "Leaf green intensity",
  leaf_base_shape        = "Leaf base shape",
  fruit_shape            = "Fruit shape",
  skin_color_intensity   = "Skin color intensity",
  flesh_color            = "Flesh color",
  ground_color           = "Ground color",
  stone_adherence        = "Stone adherence",
  shot_hole_severity     = "Shot hole severity",
  brown_rot_severity     = "Brown rot severity",
  fruit_weight           = "Fruit weight",
  petal_number           = "Petal number",
  soluble_solids_content = "Soluble solids content",
  titratable_acidity     = "Titratable acidity",
  bloom_doy              = "Full bloom date",
  dev_days               = "Fruit development period"
)
trait_label <- function(x) {
  key <- str_remove(x, "_blup$")
  out <- unname(TRAIT_LABELS[key])
  ifelse(is.na(out), x, out)
}

message("[00] 완료 | R ", getRversion(), " | PROJ_ROOT = ", PROJ_ROOT)
