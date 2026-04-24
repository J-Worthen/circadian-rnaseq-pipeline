# circadian-rnaseq-pipeline

Reproducible RNA-seq pipeline for transcript-level circadian analysis using:
SortMeRNA, STAR, RSEM, LIMBR, and ECHO.

---

## Overview

This repository contains scripts to process RNA-seq data from raw FASTQ files through transcript-level quantification and circadian analysis.

Workflow:

FASTQ → SortMeRNA → STAR → RSEM → LIMBR → ECHO

---

## Example Data

Small paired-end FASTQ files are provided for testing:

example_data/sample_R1.fastq.gz  
example_data/sample_R2.fastq.gz  

---

## Required Databases (NOT included)

Download rRNA databases from the official SortMeRNA repository:

https://github.com/sortmerna/sortmerna/tree/master/data/rRNA_databases

Required files:
- silva-euk-18s-id95.fasta
- silva-euk-28s-id98.fasta

Build indices using:

```bash
indexdb_rna \
  --ref silva-euk-18s-id95.fasta,18s_db \
  --ref silva-euk-28s-id98.fasta,28s_db
