#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Quantify transcript abundance with RSEM from paired-end,
# rRNA-removed FASTQ files.
#
# Requirements:
#   - RSEM
#   - STAR, because --star is used below
#
# Sample sheet:
#   metadata/rsem_samples.csv
#   columns: sample_id,fastq_r1,fastq_r2
# ============================================================

THREADS=8
PROJECT_DIR="$(pwd)"

FASTQ_DIR="${PROJECT_DIR}/example_data_rrna_removed"
OUTDIR="${PROJECT_DIR}/results/rsem"
REF_PREFIX="${PROJECT_DIR}/references/rsem_GRCm38_100/mouse_GRCm38_100"
SAMPLE_SHEET="${PROJECT_DIR}/metadata/rsem_samples.csv"

mkdir -p "${OUTDIR}"

[[ -f "${SAMPLE_SHEET}" ]] || { echo "ERROR: Missing sample sheet: ${SAMPLE_SHEET}" >&2; exit 1; }

# A simple check for prepared RSEM reference files.
[[ -f "${REF_PREFIX}.grp" ]] || {
  echo "ERROR: RSEM reference not found at prefix: ${REF_PREFIX}" >&2
  echo "Run scripts/04_rsem_prepare_reference.sh first or update REF_PREFIX." >&2
  exit 1
}

tail -n +2 "${SAMPLE_SHEET}" | while IFS=',' read -r SAMPLE_ID FASTQ_R1 FASTQ_R2; do
  [[ -n "${SAMPLE_ID}" ]] || continue

  R1="${FASTQ_DIR}/${FASTQ_R1}"
  R2="${FASTQ_DIR}/${FASTQ_R2}"

  for f in "${R1}" "${R2}"; do
    [[ -f "${f}" ]] || { echo "ERROR: Missing FASTQ for ${SAMPLE_ID}: ${f}" >&2; exit 1; }
  done

  echo "========================================"
  echo "Running RSEM for: ${SAMPLE_ID}"
  echo "========================================"

  rsem-calculate-expression \
    --paired-end \
    --star \
    --star-gzipped-read-file \
    -p "${THREADS}" \
    "${R1}" \
    "${R2}" \
    "${REF_PREFIX}" \
    "${OUTDIR}/${SAMPLE_ID}"

  echo "Completed: ${SAMPLE_ID}"
done

echo "All RSEM samples completed."
