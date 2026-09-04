# =============================================================================
# run_all.R — 전체 분석 파이프라인 일괄 실행
#
# 사용법
#   1) 작업 디렉터리를 프로젝트 루트("매실 유전자원")로 설정
#   2) source("R/run_all.R", encoding = "UTF-8")
#
#   또는 터미널에서:  Rscript R/run_all.R
#
# 최초 실행 시 R/00_setup.R 의 install_missing() 주석을 한 번 해제해 패키지 설치.
# =============================================================================
#install.packages("ordinal")
t0 <- Sys.time()
options(warn = 1)

STEPS <- c("00_setup.R", "01_preprocess.R", "02_diversity.R",
           "03_multivariate.R", "04_organic_acid.R")

for (s in STEPS) {
  f <- file.path("R", s)
  message("\n", strrep("=", 70), "\n>> ", s, "\n", strrep("=", 70))
  source(f, encoding = "UTF-8", local = FALSE)
}

message("\n", strrep("=", 70))
message("전체 완료 | 소요 ", round(difftime(Sys.time(), t0, units = "mins"), 1), "분")
message("표: output/tables/  |  그림: output/figures/")
message(strrep("=", 70))
