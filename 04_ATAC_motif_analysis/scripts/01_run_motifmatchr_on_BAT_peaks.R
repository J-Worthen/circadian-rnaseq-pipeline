#!/usr/bin/env Rscript

# ====================== Motif matching on BAT rhythmic peaks ======================
# Inputs:
#   - inputs/BAT_rhythmic_peaks_with_LAG_and_ZTbin.txt  (tab; has Chr, Start, End)
#   - motifs/JASPAR_circadian_TFs.meme                       (MEME v4; letter-probability PFMs)
#
# Outputs (results_motifmatchr/):
#   - motif_hits_long.tsv        (PeakID, Motif, Hit)
#   - motif_scores_long.tsv      (PeakID, Motif, Score)
#   - peaks_with_anyMotif.tsv    (input + any_motif_hit)
#   - hits_<MotifID>.bed         (per-motif BEDs of PEAKS that contain ≥1 hit)
#   - motif_hits_wide.tsv        (optional dense wide logical matrix)
#   - motif_scores_wide.tsv      (optional dense wide numeric matrix)
#   - summaries.txt              (QC summaries)
# ==================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(universalmotif)
  library(TFBSTools)
  library(motifmatchr)
  library(BSgenome.Mmusculus.UCSC.mm10)   # <--- change if not mm10
  ## NOTE:
# This script assumes mm10 genome coordinates.
# If using mm39 or another build, update the BSgenome package accordingly.
  library(Biostrings)
  library(rtracklayer)
  library(Matrix)
})

# ----------------------- user config -----------------------
project_dir <- getwd()
peaks_path <- file.path(project_dir, "inputs", "BAT_rhythmic_peaks_with_LAG_and_ZTbin.txt")
meme_path  <- file.path(project_dir, "motifs", "JASPAR_circadian_TFs.meme")
out_dir    <- file.path(project_dir, "motifmatchr_peak_hits")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# genome object (swap to mm39/hg38 if needed)
genome_obj <- BSgenome.Mmusculus.UCSC.mm10

# scanning params
p_cutoff <- 1e-03                 # tighten/relax as needed
write_per_motif_beds <- TRUE      # writes PEAK-LEVEL BEDs when TRUE
write_dense_wide      <- FALSE     # set TRUE to also write dense wide TSVs

# input coordinate system (your file looks 1-based inclusive already)
input_is_zero_based_bed <- FALSE   # set TRUE only for true BED inputs

for (f in c(peaks_path, meme_path)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

# ----------------------- read peaks ------------------------
message("Reading peaks: ", peaks_path)
dt <- data.table::fread(peaks_path, sep = "\t", header = TRUE, data.table = TRUE)

req <- c("Chr","Start","End")
if (!all(req %in% names(dt))) {
  stop("Peaks file must contain columns: ", paste(req, collapse = ", "))
}

# normalize types
dt[, Start := as.integer(Start)]
dt[, End   := as.integer(End)]

# if given BED (0-based), convert to 1-based inclusive for GRanges
if (isTRUE(input_is_zero_based_bed)) {
  dt[, Start := Start + 1L]
}

# ensure PeakID
if (!("PeakID" %in% names(dt)) || anyDuplicated(dt$PeakID)) {
  dt[, PeakID := paste0(Chr, ":", Start, "-", End)]
}

# build GRanges
message("Building GRanges…")
gr <- GRanges(
  seqnames = dt$Chr,
  ranges   = IRanges(start = dt$Start, end = dt$End),
  strand   = "*",
  PeakID   = dt$PeakID
)
stopifnot(length(gr) > 0)

# ----------------------- MEME -> PFMatrixList -----------------------
message("Reading MEME motifs: ", meme_path)
um_list <- universalmotif::read_meme(meme_path)

# keep DNA motifs only (defensive)
um_list <- Filter(function(m) m@alphabet == "DNA", um_list)
if (!length(um_list)) stop("No DNA motifs found in MEME file.")

# convert via universalmotif to TFBSTools PFMatrix objects
pf_list <- universalmotif::convert_motifs(um_list, "TFBSTools-PFMatrix")
pf_list <- Filter(function(x) methods::is(x, "PFMatrix"), pf_list)
if (!length(pf_list)) stop("No PFMatrix motifs after conversion.")

# ensure unique IDs
ids <- vapply(pf_list, TFBSTools::ID, character(1))
if (anyDuplicated(ids)) {
  ids <- make.unique(ids)
  for (i in seq_along(pf_list)) pf_list[[i]]@ID <- ids[i]   # set S4 slot directly
}

# build PFMatrixList (splat the list)
pfm <- do.call(TFBSTools::PFMatrixList, pf_list)
stopifnot(inherits(pfm, "PFMatrixList"), length(pfm) > 0)
message("Loaded ", length(pfm), " motifs.")

# ----------------------- motifmatchr run -----------------------
message("Running motifmatchr…")
mm_matches <- matchMotifs(
  pfm, gr,
  genome   = genome_obj,
  out      = "matches",
  p.cutoff = p_cutoff
)
mm_scores  <- matchMotifs(
  pfm, gr,
  genome = genome_obj,
  out    = "scores"
)

# sparse matrices
hit_mat   <- motifmatchr::motifMatches(mm_matches)   # lgCMatrix
score_mat <- motifmatchr::motifScores(mm_scores)     # dgCMatrix

# ----------------------- save LONG sparse outputs -----------------------
message("Writing sparse-long outputs…")
# hits (TRUE entries only)
hit_trip <- Matrix::summary(hit_mat)  # i, j, x
hits_long <- data.table(
  PeakID = mcols(gr)$PeakID[hit_trip$i],
  Motif  = colnames(hit_mat)[hit_trip$j],
  Hit    = as.logical(hit_trip$x)
)
fwrite(hits_long, file.path(out_dir, "motif_hits_long.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

# scores (nonzero entries only)
scr_trip <- Matrix::summary(score_mat)  # i, j, x
scores_long <- data.table(
  PeakID = mcols(gr)$PeakID[scr_trip$i],
  Motif  = colnames(score_mat)[scr_trip$j],
  Score  = scr_trip$x
)
fwrite(scores_long, file.path(out_dir, "motif_scores_long.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

# any_motif_hit merged to peaks
any_hit <- data.table(PeakID = mcols(gr)$PeakID,
                      any_motif_hit = Matrix::rowSums(hit_mat) > 0)
dt_out <- merge(dt, any_hit, by = "PeakID", all.x = TRUE)
fwrite(dt_out, file.path(out_dir, "peaks_with_anyMotif.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

# ----------------------- optional: dense wide outputs -----------------------
if (isTRUE(write_dense_wide)) {
  message("Writing dense-wide TSVs…")
  hit_dense   <- as.matrix(hit_mat)
  score_dense <- as.matrix(score_mat)
  rownames(hit_dense)   <- mcols(gr)$PeakID
  rownames(score_dense) <- mcols(gr)$PeakID
  
  fwrite(data.table(PeakID = rownames(hit_dense), as.data.table(hit_dense)),
         file.path(out_dir, "motif_hits_wide.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  fwrite(data.table(PeakID = rownames(score_dense), as.data.table(score_dense)),
         file.path(out_dir, "motif_scores_wide.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
}

# ----------------------- per-motif BEDs (PEAK-LEVEL) -----------------------
# Export each motif's *peak intervals* where ≥1 hit occurred (site-level needs different code).
if (isTRUE(write_per_motif_beds)) {
  message("Exporting per-motif BEDs (peak-level)…")
  # for each motif column, which peaks have TRUE?
  motif_ids <- colnames(hit_mat)
  for (j in seq_along(motif_ids)) {
    mot <- motif_ids[j]
    rows <- which(hit_mat[, j])
    if (!length(rows)) next
    gr_hits <- gr[rows]
    mcols(gr_hits)$name <- paste0(mot, ";", mcols(gr_hits)$PeakID)
    bed_path <- file.path(out_dir, paste0("hits_", mot, ".bed"))
    rtracklayer::export(gr_hits, bed_path, format = "BED")
  }
}

# ----------------------- summaries -----------------------
summ_file <- file.path(out_dir, "summaries.txt")
message("Writing summaries: ", summ_file)
sink(summ_file)
cat("=== motifmatchr summary ===\n")
cat("Peaks:",   nrow(dt), "\n")
cat("Motifs:",  ncol(hit_mat), "\n")
cat("p.cutoff:", p_cutoff, "\n\n")

cat("Top 15 motifs by #peaks with ≥1 hit:\n")
print(head(sort(Matrix::colSums(hit_mat), decreasing = TRUE), 15)); cat("\n")

cat("Distribution: #motifs per peak\n")
print(summary(Matrix::rowSums(hit_mat))); cat("\n")
sink()

message("Done. Outputs in: ", out_dir)

# ----------------------- OPTIONAL: exact site positions (commented) -----------------------
# If/when you need precise motif site BEDs, you can compute them with Biostrings::matchPWM.
# This is slower but does not rely on unexported motifmatchr internals.
#
# bg <- c(A=.25, C=.25, G=.25, T=.25)
# pwm_list <- lapply(as.list(pfm), function(pf)
#   TFBSTools::toPWM(pf, type="log2probratio", pseudocounts=0.8, bg=bg)
# )
# seqs <- Biostrings::getSeq(genome_obj, gr)  # DNAStringSet per peak
# min_score <- "80%"  # or a numeric threshold; adjust as needed
# dir.create(file.path(out_dir, "sites_beds"), showWarnings = FALSE)
# for (k in seq_along(pwm_list)) {
#   mot <- TFBSTools::ID(pwm_list[[k]])
#   hits_gr <- GRanges()
#   for (i in seq_along(seqs)) {
#     m1 <- Biostrings::matchPWM(pwm_list[[k]]@profileMatrix, seqs[[i]], min.score = min_score)
#     if (length(m1) > 0) {
#       # shift local positions to genomic coordinates
#       shifted <- shift(m1, start(gr[i]) - 1L)
#       tmp <- GRanges(seqnames = seqnames(gr[i]), ranges = shifted, strand = "*")
#       mcols(tmp)$name <- paste0(mot, ";", mcols(gr)$PeakID[i])
#       hits_gr <- c(hits_gr, tmp)
#     }
#   }
#   if (length(hits_gr) > 0) {
#     rtracklayer::export(hits_gr, file.path(out_dir, "sites_beds", paste0("sites_", mot, ".bed")), "BED")
#   }
# }
