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

## Step 1: Remove rRNA reads (SortMeRNA)

Script:
