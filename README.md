# Circadian RNA-seq Pipeline

This repository contains example scripts used for RNA-seq preprocessing, alignment, transcript quantification, and selected downstream motif-site processing for a circadian adipocyte transcriptomics workflow.

## Repository structure

```text
scripts/
  01_sortmerna_rrna_removal.sh
  02_star_generate_index_optional.sh
  03_star_alignment.sh
  04_rsem_prepare_reference.sh
  05_rsem_quantify_samples.sh

04_ATAC_motif_analysis/
  scripts/
    01_run_motifmatchr_on_BAT_peaks.R
    02_extract_enriched_motif_sites_by_ZT.R
    03_filter_enriched_motif_sites_by_TSS_tier.R
    04_filter_enriched_tiers_to_circadian_genes.R
    05_make_early_wave_ZT5and9_tier1_targets_list.R
  motifs/
    JASPAR_circadian_TFs.meme

metadata/
  rsem_samples.csv

references/
  Mus_musculus.GRCm38.dna.primary_assembly.fa
  Mus_musculus.GRCm38.100.gtf
```

## Notes

- Paths are written relative to the repository root using `$(pwd)` or `getwd()`.
- Reference FASTA/GTF files are expected in `references/`.
- STAR genome indexing and RSEM reference preparation require substantial memory for mammalian genomes and should be run on an HPC or high-memory workstation.
- Large generated outputs, STAR indexes, RSEM references, FASTQ files, and BAM files should not be committed to GitHub.

## Motif enrichment code availability

Motif enrichment statistics were generated using a collaborator-developed workflow that is not redistributed here. The repository includes the motif set and downstream scripts used to scan motif locations, filter enriched motif sites by genomic context, restrict sites to circadian genes, and summarize early-wave target genes from the enrichment output tables.
