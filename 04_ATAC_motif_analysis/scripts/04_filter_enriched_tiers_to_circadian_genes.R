#!/usr/bin/env Rscript

# Filter ENRICHED tier outputs to circadian genes (separate tier folders)
suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(readr)
})

project_dir  <- getwd()
in_base     <- file.path(project_dir, "motif_positions_ENRICHED_byTSS")
out_base    <- file.path(project_dir, "motif_positions_ENRICHED_byTSS_CIRCADIAN")
genes_file  <- file.path(project_dir, "files", "IVDAs_2025_ECHO_and_Metacycle_0.005_Protein_Coding_Only.txt")

stopifnot(dir.exists(in_base))
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(genes_file))

# ---- circadian gene list (col1=gene, col2=ct) ----
gdf <- suppressWarnings(read.delim(genes_file, header = TRUE, stringsAsFactors = FALSE))
if (ncol(gdf) < 2) stop("Circadian gene file must have at least 2 columns: gene, ct")
colnames(gdf)[1:2] <- c("gene", "ct")
suppressWarnings(gdf$ct <- as.numeric(gdf$ct))
circ_genes_upper <- unique(toupper(trimws(gdf$gene[!is.na(gdf$ct) & gdf$ct > 0])))
if (!length(circ_genes_upper)) stop("No circadian genes with ct > 0 found.")

# helpers
`%||%` <- function(a,b) if (!is.null(a)) a else b

infer_tier_code <- function(x) {
  # accepts "..._T1_per_ZTbin.csv" or "..._TIER1_per_ZTbin.csv"
  m <- regmatches(x, regexec("_(T(?:IER)?[1-4])_per_ZTbin\\.(csv|xlsx)$", x))[[1]]
  if (length(m) >= 2) return(sub("^TIER","T", m[2])) else return(NA_character_)  # return "T1".."T4"
}

read_tier_table <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.xlsx(path, sheet = 1)
  }
}

# find tier folders and tier files
tier_dirs <- list.dirs(in_base, full.names = TRUE, recursive = FALSE)
if (!length(tier_dirs)) stop("No tier folders found in: ", in_base)

summary_list <- list()

for (tdir in tier_dirs) {
  files_csv  <- list.files(tdir, pattern = "^motif_sites_ENRICHED_(?:T(?:IER)?[1-4])_per_ZTbin\\.csv$",
                           full.names = TRUE, ignore.case = TRUE)
  files_xlsx <- list.files(tdir, pattern = "^motif_sites_ENRICHED_(?:T(?:IER)?[1-4])_per_ZTbin\\.xlsx$",
                           full.names = TRUE, ignore.case = TRUE)
  
  tier_file <- (files_csv[1] %||% files_xlsx[1])
  if (is.null(tier_file)) {
    message("No tier CSV/XLSX in: ", tdir, " — skipping")
    next
  }
  
  tcode <- infer_tier_code(basename(tier_file))  # "T1".."T4"
  if (is.na(tcode)) {
    message("Could not infer TIER code from ", basename(tier_file), " — skipping")
    next
  }
  
  df <- read_tier_table(tier_file)
  if (!"Gene.Name" %in% names(df)) stop("Expected 'Gene.Name' in: ", tier_file)
  
  sanitize_symbol <- function(x) {
    x <- as.character(x %||% "")
    x <- sub("[|;].*$", "", x)  # keep first symbol if multiple
    toupper(trimws(x))
  }
  df$GENE_UPPER <- sanitize_symbol(df$Gene.Name)
  df_circ <- df[df$GENE_UPPER %in% circ_genes_upper, , drop = FALSE]
  
  # ensure parallel out folder
  tier_out_dir <- file.path(out_base, basename(tdir))
  dir.create(tier_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # write CSV/XLSX
  csv_out  <- file.path(tier_out_dir, paste0("motif_sites_ENRICHED_CIRCADIAN_", tcode, "_per_ZTbin.csv"))
  write.csv(df_circ, csv_out, row.names = FALSE)
  
  wb <- openxlsx::createWorkbook()
  for (bn in sort(unique(df_circ$ZT_bin))) {
    sheet_df <- df_circ[df_circ$ZT_bin == bn, , drop = FALSE]
    openxlsx::addWorksheet(wb, bn)
    openxlsx::writeData(wb, bn, sheet_df)
  }
  xlsx_out <- file.path(tier_out_dir, paste0("motif_sites_ENRICHED_CIRCADIAN_", tcode, "_per_ZTbin.xlsx"))
  openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
  
  # BEDs per ZT
  if (nrow(df_circ) > 0) {
    needed <- c("match_chr","match_start","match_end","motif","score","match_strand","ZT_bin")
    missing <- setdiff(needed, names(df_circ))
    if (length(missing)) stop("Missing columns in ", tier_file, ": ", paste(missing, collapse = ", "))
    
    for (bn in sort(unique(df_circ$ZT_bin))) {
      d2 <- df_circ[df_circ$ZT_bin == bn, needed, drop = FALSE]
      if (!nrow(d2)) next
      bed <- dplyr::select(d2, match_chr, match_start, match_end, motif, score, match_strand)
      bed$match_start <- as.integer(bed$match_start) - 1L
      colnames(bed) <- c("chrom","chromStart","chromEnd","name","score","strand")
      bed_out <- file.path(tier_out_dir, paste0("motif_sites_ENRICHED_CIRCADIAN_", bn, "_", tcode, ".bed"))
      write.table(bed, bed_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
  
  # summary rows
  if (nrow(df_circ) > 0) {
    sum_df <- df_circ %>%
      group_by(ZT_bin) %>%
      summarise(n_sites = n(),
                n_peaks = dplyr::n_distinct(peak_id),
                n_genes = dplyr::n_distinct(GENE_UPPER),
                .groups = "drop") %>%
      mutate(TIER_CODE = tcode, TIER_FOLDER = basename(tdir))
    summary_list[[length(summary_list) + 1]] <- sum_df
  }
  
  message("Wrote circadian-filtered outputs for ", basename(tdir), " → ", tier_out_dir)
}

summary_all <- dplyr::bind_rows(summary_list)
write.csv(summary_all, file.path(out_base, "ENRICHED_circadian_tier_summary_counts.csv"), row.names = FALSE)

message("Done. Circadian-filtered tier outputs in: ", out_base)
