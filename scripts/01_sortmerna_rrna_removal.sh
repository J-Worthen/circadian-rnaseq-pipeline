#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Step 1: Remove residual rRNA reads using SortMeRNA
#
# This example uses paired-end FASTQ input files and writes
# rRNA-depleted paired FASTQs for downstream alignment.
#
# Requirements:
#   - SortMeRNA
#   - Python 3
#
# Edit the variables below before running.
# =========================================================

SAMPLE="example_sample"

READ1="example_data/sample_R1.fastq.gz"
READ2="example_data/sample_R2.fastq.gz"

OUTDIR="sortmerna/${SAMPLE}"
THREADS=4

# Update this path to the local SortMeRNA rRNA database directory.
DB_DIR="sortmerna-master/data/rRNA_databases"

mkdir -p "${OUTDIR}"

MERGED="${OUTDIR}/${SAMPLE}_merged_pairs.fq"
OUT_R1="${OUTDIR}/${SAMPLE}_non_rRNA_R1.fq"
OUT_R2="${OUTDIR}/${SAMPLE}_non_rRNA_R2.fq"

for f in "${READ1}" "${READ2}"; do
  [[ -f "${f}" ]] || { echo "ERROR: Missing input FASTQ: ${f}" >&2; exit 1; }
done

for db in "${DB_DIR}/silva-euk-18s-id95.fasta" "${DB_DIR}/silva-euk-28s-id98.fasta"; do
  [[ -f "${db}" ]] || { echo "ERROR: Missing SortMeRNA database: ${db}" >&2; exit 1; }
done

echo "Merging paired FASTQ files..."

python3 - <<EOF_PY
import gzip

read1 = "${READ1}"
read2 = "${READ2}"
out = "${MERGED}"

def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

with open_fastq(read1) as r1, open_fastq(read2) as r2, open(out, "w") as merged:
    while True:
        rec1 = [r1.readline() for _ in range(4)]
        rec2 = [r2.readline() for _ in range(4)]

        if not rec1[0] and not rec2[0]:
            break
        if not rec1[0] or not rec2[0]:
            raise ValueError("R1 and R2 have different numbers of records.")

        merged.writelines(rec1)
        merged.writelines(rec2)
EOF_PY

echo "Running SortMeRNA..."

sortmerna \
  --reads "${MERGED}" \
  --ref "${DB_DIR}/silva-euk-18s-id95.fasta" \
  --ref "${DB_DIR}/silva-euk-28s-id98.fasta" \
  --paired_in \
  --fastx \
  --aligned "${OUTDIR}/${SAMPLE}_rRNA" \
  --other "${OUTDIR}/${SAMPLE}_non_rRNA" \
  -a "${THREADS}"

echo "Finding SortMeRNA non-rRNA output..."

NON_RRNA=""
if [[ -f "${OUTDIR}/${SAMPLE}_non_rRNA.fq" ]]; then
  NON_RRNA="${OUTDIR}/${SAMPLE}_non_rRNA.fq"
elif [[ -f "${OUTDIR}/${SAMPLE}_non_rRNA.fastq" ]]; then
  NON_RRNA="${OUTDIR}/${SAMPLE}_non_rRNA.fastq"
else
  echo "ERROR: Could not find SortMeRNA non-rRNA output." >&2
  ls -lh "${OUTDIR}" >&2
  exit 1
fi

echo "Unmerging non-rRNA reads back into R1/R2..."

python3 - <<EOF_PY
non_rrna = "${NON_RRNA}"
out_r1 = "${OUT_R1}"
out_r2 = "${OUT_R2}"

with open(non_rrna, "r") as inp, open(out_r1, "w") as r1, open(out_r2, "w") as r2:
    record_number = 0
    while True:
        rec = [inp.readline() for _ in range(4)]
        if not rec[0]:
            break
        if record_number % 2 == 0:
            r1.writelines(rec)
        else:
            r2.writelines(rec)
        record_number += 1

if record_number % 2 != 0:
    raise ValueError("Merged non-rRNA FASTQ has an odd number of records.")
EOF_PY

gzip -f "${OUT_R1}"
gzip -f "${OUT_R2}"

echo "Done. Non-rRNA paired FASTQs:"
echo "${OUT_R1}.gz"
echo "${OUT_R2}.gz"
