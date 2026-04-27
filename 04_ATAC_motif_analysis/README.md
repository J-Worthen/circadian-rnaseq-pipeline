# ATAC-seq Motif Analysis

This folder contains shareable downstream scripts for motif-site processing and target summarization.

## Required inputs

Place the following files in `04_ATAC_motif_analysis/inputs/` unless otherwise noted:

- `BAT_rhythmic_peaks_with_LAG_and_ZTbin.txt`
- `all_motifs_adjusted_pvalue_matrix-BAT_w_hetero.csv`
- `IVDAs_2025_ECHO_and_Metacycle_0.005_Protein_Coding_Only.txt`

Place the MEME motif file in:

04_ATAC_motif_analysis/motifs/JASPAR_circadian_TFs.meme


## Scripts

1. `01_run_motifmatchr_on_BAT_peaks.R`  
   Scans BAT rhythmic peaks for motif hits using `motifmatchr`.

2. `02_extract_enriched_motif_sites_by_ZT.R`  
   Uses the adjusted enrichment matrix to extract genomic motif positions for enriched motifs within each ZT bin.

3. `03_filter_enriched_motif_sites_by_TSS_tier.R`  
   Filters enriched motif sites into distance-to-TSS tiers.

4. `04_filter_enriched_tiers_to_circadian_genes.R`  
   Restricts tiered motif-site outputs to circadian genes.

5. `05_make_early_wave_ZT5and9_tier1_targets_list.R`  
   Generates a target-gene list for early-wave Tier 1 motif sites.

## Collaborator-developed enrichment workflow

The motif enrichment statistics themselves were generated using a collaborator-developed workflow that is not redistributed here. The scripts in this folder start from the resulting enrichment tables and public/custom motif inputs.

The file results/all_motifs_adjusted_pvalue_matrix-BAT_w_hetero.csv contains motif enrichment results used for downstream analysis. The enrichment step was performed using an external pipeline that is not distributed here. Downstream processing and analysis steps are included in this repository.