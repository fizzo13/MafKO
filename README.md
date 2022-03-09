# "Transcriptional specification of fetal to adult liver  vascular functional zonation"
Code for reproducibility and downstream analysis of the single cell RNA sequencing data available in "Transcriptional specification of fetal to adult liver  vascular functional zonation" [[add link/doi]], including: i) Single cell RNA sequencing for developmental timepoints in mouse and ii) Single cell RNA sequencing of WT and Maf KO endothelial cells

## Scripts
This repository contains the scripts used for analysis of single cell RNA seq data for "Transcriptional specification of fetal to adult liver  vascular functional zonation" [[add link/doi]], including the following files:
i) [[filename]]  Seurat base analysis of single cell RNA seq data for developmental timepoints.
ii) [[filename]] RNA velocity analysis using velocyto [[ref-link]] and scVelo [[ref-link]] for developmental timepoints to infer transitions and genes driving the transcriptional changes.
iii) [[filename]] RNA velocity analysis for either WT and Maf knockout mice, for comparison of developmental trajectories of endothelial cells upon Maf knockout.

## Available files
In addition, we provide the processed loom files and the Seurat object to facilitate downstream analysis without the need to re-process the raw data (available in GEO):
i) combined.loom - Loom hd5 file containing the combined data for all the developmental timepoints used in the publication.
ii) wt.loom - Loom hd5 file containing the processed data for WT mice.
iii) mafko.loom - Loom hd5 file containing the processed data for Maf KO mice.
iv) development.rds - .RDS file containing the processed Seurat object for all developmental timepoints
v) mafko.rds - .RDS file containing the processed Seurat object for combineed WT and Maf KO mice scRNA-seq data.

### Contact
Any questions, feel free to post on Issues or reach out to:
jmg2008@med.cornell.edu or fri2002@med.cornell.edu
