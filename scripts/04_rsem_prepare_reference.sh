#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Prepare RSEM reference for mouse GRCm38 / Ensembl release 100
#
# Requirements:
#   - RSEM
#   - STAR, because --star is used below
#
# Input files:
#   - references/Mus_musculus.GRCm38.dna.primary_assembly.fa
#   - references/Mus_musculus.GRCm38.100.gtf
#
# NOTE: This step builds an internal STAR-compatible reference and
# may require substantial memory for mammalian genomes.
# ============================================================

THREADS=8
PROJECT_DIR="$(pwd)"

FASTA="${PROJECT_DIR}/references/Mus_musculus.GRCm38.dna.primary_assembly.fa"
GTF="${PROJECT_DIR}/references/Mus_musculus.GRCm38.100.gtf"

OUTDIR="${PROJECT_DIR}/references/rsem_GRCm38_100"
REF_PREFIX="${OUTDIR}/mouse_GRCm38_100"

mkdir -p "${OUTDIR}"

for f in "${FASTA}" "${GTF}"; do
  [[ -f "${f}" ]] || { echo "ERROR: Missing required file: ${f}" >&2; exit 1; }
done

echo "Preparing RSEM reference..."
echo "FASTA: ${FASTA}"
echo "GTF:   ${GTF}"
echo "Output prefix: ${REF_PREFIX}"

rsem-prepare-reference \
  --gtf "${GTF}" \
  --star \
  -p "${THREADS}" \
  "${FASTA}" \
  "${REF_PREFIX}"

echo "RSEM reference created: ${REF_PREFIX}"
