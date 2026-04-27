#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Motif-position extraction (ENRICHED motifs only) per ZT bin using motifmatchr
# Run from the 04_ATAC_motif_analysis directory.
# Outputs to: motif_positions_ENRICHED/
#   - motif_sites_ENRICHED_per_ZTbin.csv
#   - motif_sites_ENRICHED_per_ZTbin.xlsx
#   - motif_sites_ENRICHED_<ZT>.bed  (0-based BED)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TFBSTools)
  library(universalmotif)
  library(motifmatchr)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(openxlsx)
  library(S4Vectors)
})

# --------------------------- paths & params ---------------------------
project_dir <- getwd()
peaks_path <- file.path(project_dir, "files", "BAT_rhythmic_peaks_with_LAG_and_ZTbin.txt")
pfm_rds    <- file.path(project_dir, "cache", "pfm_from_MEME.rds")
meme_path  <- file.path(project_dir, "motifs", "JASPAR_selected_73TFs_PLUS_hetero.meme")
padj_csv   <- file.path(project_dir, "files", "all_motifs_adjusted_pvalue_matrix-BAT_w_hetero.csv")

alpha_adj <- 0.01   # BH threshold for "enriched"
p_cutoff  <- 1e-3   # site-calling cutoff for positions
min_width <- 6      # ignore motif hits shorter than this
out_dir   <- file.path(project_dir, "motif_positions_ENRICHED")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genome_mm10 <- BSgenome.Mmusculus.UCSC.mm10

# --------------------------- helpers ---------------------------
make_peak_id <- function(gr) paste0(seqnames(gr), ":", start(gr), "-", end(gr))

score_col <- function(gr) {
  if (length(gr) == 0) return(numeric(0))
  nm <- names(mcols(gr))
  cand <- c("score","motif.score","log10.pvalue","p.value")
  hit <- intersect(cand, nm)
  if (length(hit) == 0) return(rep(NA_real_, length(gr)))
  as.numeric(mcols(gr)[[hit[1]]])
}

# ------------------------ load peaks (GRanges) ------------------------
stopifnot(file.exists(peaks_path))
peaks_df <- read.delim(peaks_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!("Strand" %in% names(peaks_df))) peaks_df$Strand <- "*"
peaks_df$Strand[!peaks_df$Strand %in% c("+","-","*")] <- "*"

peaks_gr <- GRanges(
  seqnames = Rle(peaks_df$Chr),
  ranges   = IRanges(start = as.integer(peaks_df$Start), end = as.integer(peaks_df$End)),
  strand   = Rle(peaks_df$Strand)
)
mcols(peaks_gr) <- peaks_df[, setdiff(names(peaks_df), c("Chr","Start","End","Strand")), drop = FALSE]
stopifnot("ZT_bin" %in% names(mcols(peaks_gr)))
peaks_gr_l <- split(peaks_gr, mcols(peaks_gr)$ZT_bin)

# ------------------------ load PFMs (cache or rebuild) ------------------------
if (!file.exists(pfm_rds)) {
  message("PFM cache missing at ", pfm_rds, " — rebuilding from MEME.")
  stopifnot(file.exists(meme_path))
  um <- universalmotif::read_meme(meme_path)
  um <- Filter(function(m) m@alphabet == "DNA", um)
  stopifnot(length(um) > 0)
  pf_list <- universalmotif::convert_motifs(um, "TFBSTools-PFMatrix")
  pf_list <- Filter(function(x) methods::is(x, "PFMatrix"), pf_list)
  stopifnot(length(pf_list) > 0)
  ids <- vapply(pf_list, TFBSTools::ID, character(1))
  if (anyDuplicated(ids)) {
    ids <- make.unique(ids)
    for (i in seq_along(pf_list)) pf_list[[i]]@ID <- ids[i]
  }
  pfm <- do.call(TFBSTools::PFMatrixList, pf_list)
  dir.create(dirname(pfm_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(pfm, pfm_rds)
  message("Saved PFMs to: ", pfm_rds)
}

PFMatrixList_sub <- readRDS(pfm_rds)
stopifnot(inherits(PFMatrixList_sub, "PFMatrixList"), length(PFMatrixList_sub) > 0)
names(PFMatrixList_sub) <- TFBSTools::ID(PFMatrixList_sub)

# motif -> TF mapping (for annotation)
motif_to_gene_df <- data.frame(
  motif = TFBSTools::ID(PFMatrixList_sub),
  gene  = TFBSTools::name(PFMatrixList_sub),
  stringsAsFactors = FALSE
)

# ---------------- load enrichment padj and pick enriched motifs ----------------
stopifnot(file.exists(padj_csv))
padj_df <- read.csv(padj_csv, check.names = FALSE, row.names = 1)

# rownames like "MA0603.2::Arntl" → extract motif ID before '::'
motif_ids_from_rows <- sub("::.*$", "", rownames(padj_df))

# per-ZT enriched motif IDs
enriched_by_bin <- lapply(colnames(padj_df), function(col) {
  ok <- is.finite(padj_df[[col]]) & (padj_df[[col]] <= alpha_adj)
  unique(motif_ids_from_rows[ok])
})
names(enriched_by_bin) <- colnames(padj_df)

# intersect ZT labels between enrichment and peaks
common_bins <- intersect(names(enriched_by_bin), names(peaks_gr_l))
if (length(common_bins) == 0) stop("ZT bin names don't overlap between enrichment and peaks.")

# ---------------- core extractor for one ZT bin -------------------
extract_motif_positions_for_bin <- function(peaks_bin_gr, bin_name,
                                            pwms, genome, p.cutoff = 1e-3, min_width = 6) {
  if (length(peaks_bin_gr) == 0) return(data.frame())
  peak_id <- make_peak_id(peaks_bin_gr)
  
  pos_by_motif <- motifmatchr::matchMotifs(
    pwms, peaks_bin_gr, genome = genome, out = "positions", p.cutoff = p.cutoff
  )
  
  motif_hits_list <- lapply(names(pos_by_motif), function(mname) {
    pos_obj <- pos_by_motif[[mname]]
    if (length(pos_obj) == 0) return(NULL)
    
    if (inherits(pos_obj, "GRangesList")) {
      as_df <- function(gr, i) {
        if (length(gr) == 0) return(NULL)
        data.frame(
          peak_ix       = i,
          match_chr     = as.character(seqnames(gr)),
          match_start   = start(gr),
          match_end     = end(gr),
          match_strand  = as.character(strand(gr)),
          score         = score_col(gr),
          stringsAsFactors = FALSE
        )
      }
      tmp <- Map(as_df, as.list(pos_obj), seq_along(pos_obj))
      res <- if (length(tmp)) do.call(rbind, tmp) else NULL
      
    } else if (inherits(pos_obj, "GRanges")) {
      ol <- GenomicRanges::findOverlaps(pos_obj, peaks_bin_gr, ignore.strand = TRUE)
      if (length(ol) == 0) return(NULL)
      qh <- S4Vectors::queryHits(ol)
      sh <- S4Vectors::subjectHits(ol)
      sc <- score_col(pos_obj)
      res <- data.frame(
        peak_ix       = sh,
        match_chr     = as.character(seqnames(pos_obj))[qh],
        match_start   = start(pos_obj)[qh],
        match_end     = end(pos_obj)[qh],
        match_strand  = as.character(strand(pos_obj))[qh],
        score         = if (length(sc)) sc[qh] else NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      return(NULL)
    }
    
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$motif <- mname
    res
  })
  
  hits <- if (length(motif_hits_list)) do.call(rbind, motif_hits_list) else NULL
  if (is.null(hits) || nrow(hits) == 0) return(data.frame())
  
  # motif -> TF
  hits <- merge(hits, motif_to_gene_df, by = "motif", all.x = TRUE)
  
  # parent peak info
  peak_df <- data.frame(
    peak_ix      = seq_along(peaks_bin_gr),
    peak_chr     = as.character(seqnames(peaks_bin_gr)),
    peak_start   = start(peaks_bin_gr),
    peak_end     = end(peaks_bin_gr),
    peak_strand  = as.character(strand(peaks_bin_gr)),
    peak_id      = peak_id,
    stringsAsFactors = FALSE
  )
  if (ncol(mcols(peaks_bin_gr)) > 0) {
    meta_df <- as.data.frame(mcols(peaks_bin_gr))
    colnames(meta_df) <- make.unique(colnames(meta_df))
    peak_df <- cbind(peak_df, meta_df)
  }
  
  hits2 <- merge(hits, peak_df, by = "peak_ix", all.x = TRUE)
  hits2 <- hits2[(hits2$match_end - hits2$match_start + 1) >= min_width, , drop = FALSE]
  if (!nrow(hits2)) return(data.frame())
  hits2$ZT_bin <- bin_name
  
  keep_cols <- c(
    "ZT_bin", "motif", "gene",
    "peak_id", "peak_chr", "peak_start", "peak_end", "peak_strand",
    "match_chr", "match_start", "match_end", "match_strand",
    "score"
  )
  extra_cols <- setdiff(colnames(hits2), c(keep_cols, "peak_ix"))
  hits2 <- hits2[, c(keep_cols, extra_cols), drop = FALSE]
  hits2[order(hits2$motif, hits2$peak_chr, hits2$match_start), ]
}

# ---------------------- run bins (ENRICHED motifs only) ----------------------
message("Scanning ZT bins (ENRICHED motifs only; BH ≤ ", alpha_adj,
        ", p_cutoff = ", p_cutoff, ", min_width = ", min_width, ")")

motif_hits_by_bin <- lapply(common_bins, function(bn) {
  grb <- peaks_gr_l[[bn]]
  if (length(grb) == 0) {
    message("  - ", bn, ": 0 peaks → skip.")
    return(data.frame())
  }
  
  mot_ids <- enriched_by_bin[[bn]]
  if (is.null(mot_ids) || length(mot_ids) == 0) {
    message("  - ", bn, ": no enriched motifs → skip.")
    return(data.frame())
  }
  
  mot_ids <- intersect(mot_ids, names(PFMatrixList_sub))
  if (length(mot_ids) == 0) {
    message("  - ", bn, ": enriched motif IDs not present in PFMs → skip.")
    return(data.frame())
  }
  
  pwms_sub <- PFMatrixList_sub[mot_ids]
  message("  - ", bn, ": ", length(grb), " peaks; ", length(pwms_sub), " enriched motifs")
  extract_motif_positions_for_bin(
    peaks_bin_gr = grb,
    bin_name     = bn,
    pwms         = pwms_sub,
    genome       = genome_mm10,
    p.cutoff     = p_cutoff,
    min_width    = min_width
  )
})
names(motif_hits_by_bin) <- common_bins

motif_hits_all <- do.call(rbind, motif_hits_by_bin)
message("Total ENRICHED motif matches kept: ", ifelse(is.null(nrow(motif_hits_all)), 0L, nrow(motif_hits_all)))

# ---------------------- write CSV + XLSX (ENRICHED) ----------------------
csv_out  <- file.path(out_dir, "motif_sites_ENRICHED_per_ZTbin.csv")
xlsx_out <- file.path(out_dir, "motif_sites_ENRICHED_per_ZTbin.xlsx")

if (!is.null(motif_hits_all) && nrow(motif_hits_all) > 0) {
  write.csv(motif_hits_all, csv_out, row.names = FALSE)
  wb_sites <- openxlsx::createWorkbook()
  for (bn in names(motif_hits_by_bin)) {
    openxlsx::addWorksheet(wb_sites, bn)
    openxlsx::writeData(wb_sites, bn, motif_hits_by_bin[[bn]])
  }
  openxlsx::saveWorkbook(wb_sites, xlsx_out, overwrite = TRUE)
} else {
  message("No ENRICHED motif sites found at the current thresholds; skipping CSV/XLSX write.")
}

# ---------------------- export 0-based BED per bin (ENRICHED) ----------------------
for (bn in names(motif_hits_by_bin)) {
  df <- motif_hits_by_bin[[bn]]
  if (!is.null(df) && nrow(df) > 0) {
    bed <- df[, c("match_chr", "match_start", "match_end", "motif", "score", "match_strand")]
    bed$match_start <- as.integer(bed$match_start) - 1L  # 0-based BED
    colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
    fn <- file.path(out_dir, paste0("motif_sites_ENRICHED_", bn, ".bed"))
    write.table(bed, fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

message("Done. ENRICHED outputs in: ", out_dir)
