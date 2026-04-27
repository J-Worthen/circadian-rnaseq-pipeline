#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Step 2: Align rRNA-removed paired-end RNA-seq reads with STAR
#
# This script assumes that a valid STAR genome index already exists.
# Generate the index separately with 02_star_generate_index_optional.sh
# or provide a compatible prebuilt STAR index.
# =========================================================

SAMPLE="example_sample"

READ1="example_data_rrna_removed/sample_clean_R1.fastq.gz"
READ2="example_data_rrna_removed/sample_clean_R2.fastq.gz"

PROJECT_DIR="$(pwd)"
STAR_INDEX="${PROJECT_DIR}/references/star_index"
OUTDIR="${PROJECT_DIR}/star/${SAMPLE}"
THREADS=8

mkdir -p "${OUTDIR}"

echo "Checking STAR installation..."
STAR --version

echo "Checking input FASTQs..."
for f in "${READ1}" "${READ2}"; do
  [[ -f "${f}" ]] || { echo "ERROR: Missing input FASTQ: ${f}" >&2; exit 1; }
done

echo "Checking STAR index..."
for f in "Genome" "SA" "SAindex"; do
  [[ -f "${STAR_INDEX}/${f}" ]] || {
    echo "ERROR: STAR index appears incomplete. Missing: ${STAR_INDEX}/${f}" >&2
    echo "Generate the index first or update STAR_INDEX to a valid prebuilt index." >&2
    exit 1
  }
done

echo "Running STAR alignment for ${SAMPLE}..."

STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn "${READ1}" "${READ2}" \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --twopassMode Basic \
  --outFilterMultimapNmax 1 \
  --quantMode TranscriptomeSAM \
  --limitBAMsortRAM 16000000000 \
  --outFileNamePrefix "${OUTDIR}/${SAMPLE}_"

echo "Done. Coordinate-sorted BAM:"
echo "${OUTDIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
echo "Transcriptome BAM for transcript quantification:"
echo "${OUTDIR}/${SAMPLE}_Aligned.toTranscriptome.out.bam"
