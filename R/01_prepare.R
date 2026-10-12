# =============================================================================
# 01_prepare.R  Curation, exclusions, and stage-1 within-year adjustment
# =============================================================================
raw <- read_csv(file.path(ROOT, "data", "phenotypes.csv"), col_types = cols(.default = col_character()),
                na = c("", "NA"), show_col_types = FALSE)
names(raw)[1] <- sub("^﻿", "", names(raw)[1])
stopifnot(!anyDuplicated(paste(raw$no, raw$year)))
raw[] <- lapply(raw, function(x) ifelse(trimws(x) == "", NA_character_, trimws(x)))
dat <- as.data.frame(raw)
dat$year <- as.integer(dat$year)

# ---- exclusions -------------------------------------------------------------
n_before <- nrow(dat)
dat <- dat[!dat$no %in% EXCLUDE_ACCESSIONS, ]
dat <- dat[, setdiff(names(dat), EXCLUDE_TRAITS)]
message("Excluded ", n_before - nrow(dat), " records of ", paste(EXCLUDE_ACCESSIONS, collapse = ", "),
        " and the descriptors ", paste(EXCLUDE_TRAITS, collapse = ", "))
# varietal characteristics: one value per accession in every year of record
for (tr in FIXED_TRAITS) stopifnot(all(tapply(dat[[tr]], dat$no, function(v) length(unique(na.omit(v))) <= 1)))

# ---- curation ---------------------------------------------------------------
for (tr in c("fruit_weight", "petal_number", "soluble_solids_content", "titratable_acidity")) dat[[tr]] <- as.numeric(dat[[tr]])
dat$fruit_weight[dat$year == 2024] <- dat$fruit_weight[dat$year == 2024] / 3    # 2024 values were three-fruit totals
dat$skin_color_intensity <- recode(dat$skin_color_intensity, stong = "strong", `very stong` = "very strong")
dat$leaf_green_intensity[dat$leaf_green_intensity == "other"] <- NA_character_
for (tr in ORD) stopifnot(all(is.na(dat[[tr]]) | dat[[tr]] %in% LEV[[tr]]))

bloom_date   <- as.Date(dat$full_bloom_date); harvest_date <- as.Date(dat$maturity_date)
invalid_bloom   <- !is.na(bloom_date)   & (format(bloom_date, "%m") > "05"   | as.integer(format(bloom_date, "%Y"))   != dat$year)
invalid_harvest <- !is.na(harvest_date) & (format(harvest_date, "%m") > "08" | as.integer(format(harvest_date, "%Y")) != dat$year)
date_log <- bind_rows(
  data.frame(no = dat$no[invalid_bloom],   year = dat$year[invalid_bloom],   field = rep("full_bloom_date", sum(invalid_bloom)),   original = dat$full_bloom_date[invalid_bloom]),
  data.frame(no = dat$no[invalid_harvest], year = dat$year[invalid_harvest], field = rep("harvest_date",    sum(invalid_harvest)), original = dat$maturity_date[invalid_harvest]))
save_table(date_log, "date_curation_log")
dat$bloom_doy_recorded <- as.integer(format(bloom_date, "%j"))      # recorded value, kept in the curated records for reference
bloom_date[invalid_bloom] <- NA; harvest_date[invalid_harvest] <- NA
dat$bloom_doy   <- as.integer(format(bloom_date, "%j"))
dat$harvest_doy <- as.integer(format(harvest_date, "%j"))
dat$batch <- ifelse(is.na(harvest_date), NA_character_, paste(dat$year, as.character(harvest_date), sep = "_"))
dat$no <- as.character(dat$no)

# ---- stage 1: within-year adjustment ---------------------------------------
# Fruit traits: subtract the median of the harvest batch (year x date) and add the overall median back,
# so that adjusted values keep the original unit. Other quantitative traits: the same with the year median.
grand_med <- sapply(CONT, function(tr) median(dat[[tr]], na.rm = TRUE))
adjust <- function(y, group, gm) { med <- ave(y, group, FUN = function(v) median(v, na.rm = TRUE)); y - med + gm }
for (tr in CONT) {
  g <- if (tr %in% FRUIT_CONT) dat$batch else as.character(dat$year)
  dat[[paste0(tr, "_adj")]] <- adjust(dat[[tr]], g, grand_med[[tr]])
}
save_table(dat, "curated_adjusted_records")

# ---- descriptive tables -----------------------------------------------------
years_by_id <- dat %>% count(no, name = "n_years") %>% arrange(no)
save_table(years_by_id, "observation_years_by_accession")
save_table(dat %>% count(year, name = "n_records"), "records_by_year")
save_table(bind_rows(lapply(PRIMARY, function(tr) dat %>% group_by(year) %>%
  summarise(n = sum(!is.na(.data[[tr]])), .groups = "drop") %>% mutate(trait = tr))), "observations_by_trait_year")
save_table(dat %>% select(no, year, all_of(CONT)) %>% pivot_longer(-c(no, year), names_to = "trait", values_to = "value") %>%
  group_by(trait, year) %>% summarise(n = sum(is.finite(value)), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
                                      minimum = min(value, na.rm = TRUE), maximum = max(value, na.rm = TRUE), .groups = "drop"),
  "annual_continuous_summary")
# harvest batches with their medians (shows the size of the stage-1 adjustment)
save_table(dat %>% filter(!is.na(batch)) %>% group_by(year, batch, harvest_date = as.character(harvest_date)) %>%
  summarise(n_accessions = n(), fruit_weight_median = median(fruit_weight, na.rm = TRUE),
            ssc_median = median(soluble_solids_content, na.rm = TRUE), ta_median = median(titratable_acidity, na.rm = TRUE),
            .groups = "drop"), "harvest_batches")
save_table(dat %>% select(no, year, all_of(CAT)) %>% pivot_longer(-c(no, year), names_to = "trait", values_to = "level") %>%
  filter(!is.na(level)) %>% count(trait, level, name = "n_records"), "category_frequencies")
message("Prepared ", nrow(dat), " accession-year records, ", n_distinct(dat$no), " accessions, ", length(PRIMARY), " primary variables")
