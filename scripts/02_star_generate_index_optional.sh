#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Optional: Generate a STAR genome index
#
# NOTE: STAR indexing for mammalian genomes requires substantial
# memory and should be run on an HPC or high-memory workstation.
# This step may be slow or fail on standard laptops.
# =========================================================

THREADS=8
PROJECT_DIR="$(pwd)"

GENOME_FASTA="${PROJECT_DIR}/references/Mus_musculus.GRCm38.dna.primary_assembly.fa"
GTF="${PROJECT_DIR}/references/Mus_musculus.GRCm38.100.gtf"
STAR_INDEX="${PROJECT_DIR}/references/star_index"

# For 150 bp paired-end reads, use 149. For 100 bp reads, use 99.
SJDB_OVERHANG=149

mkdir -p "${STAR_INDEX}"

for f in "${GENOME_FASTA}" "${GTF}"; do
  [[ -f "${f}" ]] || { echo "ERROR: Missing required file: ${f}" >&2; exit 1; }
done

STAR --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang "${SJDB_OVERHANG}"

echo "STAR index written to: ${STAR_INDEX}"
