#!/usr/bin/env Rscript

# Create a gene list for early-wave Tier 1 BAT motif targets (ZT5/ZT9)
# that overlap IVDA circadian genes.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(purrr)
})

# ===================== CONFIG =====================
tier1_csv   <- file.path("motif_positions_ENRICHED_byTSS_CIRCADIAN", "TIER1_promoter__2kb", "motif_sites_ENRICHED_CIRCADIAN_T1_per_ZTbin.csv")
ivda_ct_txt <- file.path("files", "IVDAs_2025_ECHO_and_Metacycle_0.005_Protein_Coding_Only.txt")  # 2 cols, no header
out_dir     <- "files"
out_csv     <- file.path(out_dir, "circadian_early_wave_targets.csv")

# ===================== HELPERS =====================
# case-insensitive column picker; returns the first matching name or NA
pick_col <- function(df, pattern_vec) {
  nms <- names(df)
  m <- map_lgl(nms, ~ any(str_detect(tolower(.x), tolower(pattern_vec))))
  if (!any(m)) return(NA_character_)
  nms[which(m)[1]]
}

# Direct name match (regex), case-insensitive; returns vector of matching names
grep_cols <- function(df, regex) {
  nms <- names(df)
  nms[str_detect(tolower(nms), tolower(regex))]
}

norm_sym <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x <- gsub("[[:space:]\u200B\u200C\u200D]+", "", x, perl = TRUE) # strip weird spaces
  x
}

# ===================== LOAD TIER1 =====================
if (!file.exists(tier1_csv)) stop("Missing Tier-1 CSV: ", tier1_csv)
t1 <- suppressMessages(readr::read_csv(tier1_csv, show_col_types = FALSE))
cat("Tier-1 rows:", nrow(t1), " | columns:", length(names(t1)), "\n")

# ===================== FILTER: TIER 1 ONLY =====================
col_tiercode <- grep_cols(t1, "^TIER[_ ]?CODE$|^tier[_ ]?code$")
col_tsstier  <- grep_cols(t1, "^TSS[_ ]?TIER$|^tss[_ ]?tier$")

is_tier1 <- rep(TRUE, nrow(t1))
if (length(col_tiercode)) {
  is_tier1 <- is_tier1 & toupper(t1[[col_tiercode[1]]]) == "T1"
}
if (length(col_tsstier)) {
  is_tier1 <- is_tier1 & str_detect(toupper(t1[[col_tsstier[1]]]), "TIER1")
}

t1 <- t1[is_tier1, , drop = FALSE]
cat("After Tier-1 filter:", nrow(t1), "rows\n")

# ===================== EARLY-WAVE: ZT5/ZT9 replicate columns (>0) =====================
zt5_cols <- grep_cols(t1, "^ZT\\s*0?5_.*_BAT$")
zt9_cols <- grep_cols(t1, "^ZT\\s*0?9_.*_BAT$")

cat("Detected ZT5 columns:", paste(zt5_cols, collapse = ", "), "\n")
cat("Detected ZT9 columns:", paste(zt9_cols, collapse = ", "), "\n")

if (length(zt5_cols) == 0 && length(zt9_cols) == 0) {
  stop("No early-wave columns found (expected names like ZT5_1_BAT, ZT9_2_BAT).")
}

numify <- function(x) suppressWarnings(as.numeric(as.character(x)))
early_mask <- rep(FALSE, nrow(t1))
if (length(zt5_cols)) {
  early_mask <- early_mask | rowSums(as.data.frame(lapply(t1[zt5_cols], numify)) %>% replace(is.na(.), 0)) > 0
}
if (length(zt9_cols)) {
  early_mask <- early_mask | rowSums(as.data.frame(lapply(t1[zt9_cols], numify)) %>% replace(is.na(.), 0)) > 0
}

t1_early <- t1[early_mask, , drop = FALSE]
cat("After early-wave (ZT5/9) filter:", nrow(t1_early), "rows\n")

# ===================== GENE SYMBOL (prefer GENE_UPPER) =====================
col_gene_upper <- grep_cols(t1_early, "^GENE[_ ]?UPPER$")
col_gene_fallback <- c(
  grep_cols(t1_early, "^Gene[.]?Name(\\.{3}\\d+)?$"),
  grep_cols(t1_early, "^gene$"),
  grep_cols(t1_early, "^gene[_ ]?symbol$")
)

if (length(col_gene_upper)) {
  genes <- norm_sym(t1_early[[col_gene_upper[1]]])
} else if (length(col_gene_fallback)) {
  genes <- norm_sym(t1_early[[col_gene_fallback[1]]])
} else {
  stop("Could not find a gene symbol column (GENE_UPPER / Gene.Name / gene / gene_symbol).")
}

early_genes <- tibble(gene_symbol = genes) %>%
  filter(nzchar(gene_symbol)) %>%
  distinct()

cat("Unique early-wave Tier-1 genes (pre-CT join):", nrow(early_genes), "\n")

# ===================== LOAD IVDA CTs =====================
if (!file.exists(ivda_ct_txt)) stop("Missing IVDA CT file: ", ivda_ct_txt)
ivda <- suppressMessages(readr::read_tsv(ivda_ct_txt, col_names = FALSE, show_col_types = FALSE))
if (ncol(ivda) < 2) stop("IVDA CT file must have 2 columns: gene_symbol, CT (no header).")

ivda <- ivda %>%
  transmute(
    gene_symbol = norm_sym(X1),
    ct          = suppressWarnings(as.numeric(X2))
  ) %>%
  filter(nzchar(gene_symbol), !is.na(ct)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

cat("IVDA genes with CT:", nrow(ivda), "\n")

# ===================== INTERSECT & WRITE =====================
out_tbl <- early_genes %>%
  inner_join(ivda, by = "gene_symbol") %>%
  arrange(gene_symbol)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
readr::write_csv(out_tbl %>% select(gene_symbol, ct), out_csv)

cat("Wrote early-wave list:\n - ", out_csv, "\n", sep = "")
cat("Summary: ", nrow(out_tbl), " genes (unique Tier1 early-wave targets with CT)\n", sep = "")
