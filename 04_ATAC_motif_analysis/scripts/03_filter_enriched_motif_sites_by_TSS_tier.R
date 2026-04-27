#!/usr/bin/env Rscript
# ==============================================================================
# Post-filter enriched motif sites by distance-to-TSS tiers (separate folders)
# - Tiers by |Distance.to.TSS|:
#     T1: ≤ 2 kb
#     T2: 2–10 kb
#     T3: 10–100 kb
#     T4: > 100 kb
# - Constraint: T2–T4 must be UPSTREAM (Distance.to.TSS < 0) only
#
# INPUT  : motif_positions_ENRICHED*/motif_sites_ENRICHED*_per_ZTbin.csv
# OUTPUT : motif_positions_ENRICHED*_byTSS/
#   ├─ TIER1_promoter__2kb/
#   │    ├─ motif_sites_ENRICHED_T1_per_ZTbin.csv
#   │    ├─ motif_sites_ENRICHED_T1_per_ZTbin.xlsx
#   │    ├─ motif_sites_ENRICHED_<ZT>_T1.bed
#   ├─ TIER2_proximal_2to10kb/   (UPSTREAM only)
#   ├─ TIER3_distal_10to100kb/   (UPSTREAM only)
#   └─ TIER4_far__100kb/         (UPSTREAM only)
# Also writes ENRICHED_tier_summary_counts.csv with n_sites/n_peaks/n_genes.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(openxlsx)
  library(janitor)
  library(stringr)
})

# ------------------------ USER PATHS (edit) ------------------------
# Example (non-circadian):
# in_csv   <- file.path(getwd(), "motif_positions_ENRICHED", "motif_sites_ENRICHED_per_ZTbin.csv")
# out_base <- file.path(getwd(), "motif_positions_ENRICHED_byTSS")

# Example (circadian):
# in_csv   <- file.path(getwd(), "motif_positions_ENRICHED_CIRCADIAN", "motif_sites_ENRICHED_CIRCADIAN_per_ZTbin.csv")
# out_base <- file.path(getwd(), "motif_positions_ENRICHED_byTSS_CIRCADIAN")

# ---- EDIT THESE TWO LINES FOR EACH RUN ----
in_csv   <- file.path(getwd(), "motif_positions_ENRICHED", "motif_sites_ENRICHED_per_ZTbin.csv")
out_base <- file.path(getwd(), "motif_positions_ENRICHED_byTSS")
# ---------------------------------------------------------------

stopifnot(file.exists(in_csv))
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# Consistent folder names (match prior convention)
tier_dirname <- c(
  T1 = "TIER1_promoter__2kb",
  T2 = "TIER2_proximal_2to10kb",
  T3 = "TIER3_distal_10to100kb",
  T4 = "TIER4_far__100kb"
)

# ---------- helpers ----------
derive_tier_code <- function(dist_bp) {
  if (is.na(dist_bp)) return(NA_character_)
  dkb <- abs(dist_bp) / 1000
  if (dkb <= 2)       return("T1")
  if (dkb <= 10)      return("T2")
  if (dkb <= 100)     return("T3")
  return("T4")
}

# prefer robust column names regardless of casing / punctuation
pick_col <- function(df, candidates) {
  got <- candidates[candidates %in% names(df)]
  if (length(got)) got[1] else NA_character_
}

# ---------- load ----------
raw <- read_csv(in_csv, show_col_types = FALSE) %>% clean_names()

# required columns
stopifnot(all(c("zt_bin","distance_to_tss") %in% names(raw)))

# optional variant columns
col_peak  <- pick_col(raw, c("peak_id","peakid"))   # peak identity
col_gene  <- pick_col(raw, c("gene_name","gene","gene_symbol","target_gene","gene_name")) # gene symbol
if (is.na(col_gene)) stop("Could not locate a gene name column in input.")

df <- raw %>%
  mutate(
    zt_bin = toupper(if_else(str_detect(zt_bin, "^ZT", negate = FALSE), zt_bin, paste0("ZT", zt_bin))),
    dist   = as.numeric(distance_to_tss),
    TIER_CODE = vapply(dist, derive_tier_code, character(1))
  ) %>%
  filter(!is.na(TIER_CODE), !is.na(zt_bin)) %>%
  # Enforce upstream-only for T2–T4
  filter(!(TIER_CODE %in% c("T2","T3","T4") & dist >= 0))

# Assign human-readable tier labels (matching folder names)
df <- df %>%
  mutate(
    TSS_TIER = dplyr::recode(TIER_CODE,
                             T1 = "TIER1_promoter__2kb",
                             T2 = "TIER2_proximal_2to10kb",
                             T3 = "TIER3_distal_10to100kb",
                             T4 = "TIER4_far__100kb",
                             .default = NA_character_
    )
  ) %>%
  filter(!is.na(TSS_TIER))

# ---------- write per-tier outputs ----------
tier_levels <- c("T1","T2","T3","T4")
zt_levels   <- c("ZT1","ZT5","ZT9","ZT13","ZT17","ZT21")

for (tcode in tier_levels) {
  tlabel <- dplyr::recode(tcode,
                          T1="TIER1_promoter__2kb",
                          T2="TIER2_proximal_2to10kb",
                          T3="TIER3_distal_10to100kb",
                          T4="TIER4_far__100kb")
  sub_df <- df %>% filter(TIER_CODE == tcode)
  
  tier_dir <- file.path(out_base, tlabel)
  dir.create(tier_dir, showWarnings = FALSE, recursive = TRUE)
  
  # CSV (tier-only)
  csv_out <- file.path(tier_dir, paste0("motif_sites_ENRICHED_", tcode, "_per_ZTbin.csv"))
  write_csv(sub_df, csv_out)
  
  # Excel workbook: one sheet per ZT bin
  wb <- openxlsx::createWorkbook()
  zts <- sort(unique(sub_df$zt_bin))
  if (length(zts) == 0) {
    # still create an empty workbook to be explicit
    openxlsx::addWorksheet(wb, "EMPTY")
    openxlsx::writeData(wb, "EMPTY", tibble())
  } else {
    for (bn in zts) {
      sheet_df <- sub_df %>% filter(zt_bin == bn)
      openxlsx::addWorksheet(wb, bn)
      openxlsx::writeData(wb, bn, sheet_df)
    }
  }
  xlsx_out <- file.path(tier_dir, paste0("motif_sites_ENRICHED_", tcode, "_per_ZTbin.xlsx"))
  openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
  
  # Per-ZT BEDs (0-based)
  if (nrow(sub_df) > 0) {
    req_cols <- c("match_chr","match_start","match_end","motif","score","match_strand")
    missing  <- setdiff(req_cols, names(sub_df))
    if (length(missing)) {
      warning("Missing columns for BED in ", tcode, ": ", paste(missing, collapse = ", "),
              ". Skipping BED writes for this tier.")
    } else {
      for (bn in sort(unique(sub_df$zt_bin))) {
        d2 <- sub_df[sub_df$zt_bin == bn, , drop = FALSE]
        if (nrow(d2) == 0) next
        bed <- d2[, req_cols]
        bed$match_start <- as.integer(bed$match_start) - 1L  # 0-based start
        colnames(bed) <- c("chrom","chromStart","chromEnd","name","score","strand")
        bed_out <- file.path(tier_dir, paste0("motif_sites_ENRICHED_", bn, "_", tcode, ".bed"))
        write.table(bed, bed_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}

# ---------- summary counts (n_sites, n_peaks, n_genes) ----------
summary_df <- df %>%
  mutate(
    peak_id_sane = if (!is.na(col_peak)) .data[[col_peak]] else NA_character_,
    gene_sym     = .data[[col_gene]]
  ) %>%
  group_by(TIER_CODE, TSS_TIER, zt_bin) %>%
  summarise(
    n_sites = n(),
    n_peaks = if (!is.na(col_peak)) n_distinct(peak_id_sane) else NA_integer_,
    n_genes = n_distinct(gene_sym),
    .groups = "drop"
  ) %>%
  arrange(TIER_CODE, zt_bin)

write_csv(summary_df, file.path(out_base, "ENRICHED_tier_summary_counts.csv"))

message("Done. Tiered outputs written under: ", normalizePath(out_base))
