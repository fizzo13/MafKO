# "Transcriptional specification of fetal to adult liver  vascular functional zonation"
Code for reproducibility and downstream analysis of the single cell RNA sequencing data available in "Transcriptional specification of fetal to adult liver  vascular functional zonation" [[add link/doi]], including: i) Single cell RNA sequencing for developmental timepoints in mouse and ii) Single cell RNA sequencing of WT and Maf KO endothelial cells

## Scripts
This repository contains the scripts used for analysis of single cell RNA seq data for "Transcriptional specification of fetal to adult liver  vascular functional zonation" [[add link/doi]], including the following files:
</p>
i) <strong>development.R</strong> - Seurat base analysis of single cell RNA seq data for developmental timepoints.
</p>
ii) <strong>scVelo_development.ipynb</strong> and <strong>scVelo_development.html</strong> RNA velocity analysis using velocyto (https://github.com/velocyto-team/velocyto.py) and scVelo (https://github.com/theislab/scvelo) for developmental timepoints to infer transitions and genes driving the transcriptional changes. The .ipynb file contains the Jupyter notebook and the .html file includes graphical output.
</p>
iii) <strong>MafKO.R</strong> Seurat base analysis for single cell RNA seq data for WT vs Maf KO experiment
</p>
iv) <strong>scVelo_MafKO.ipynb</strong> and <strong>scVelo_MafKO.html</strong> RNA velocity analysis for either WT and Maf knockout mice, for comparison of developmental trajectories of endothelial cells upon Maf knockout.The .ipynb file contains the Jupyter notebook and the .html file includes graphical output.

## Available files upon request
Raw data for this publication is available at GEO (GSE174209 and GSE174208). In addition, we can provide the processed loom files and the Seurat object to facilitate downstream analysis without the need to re-process the raw data. Due to file size constrains, we are not able to upload the files here, but the following files are available upon request:
</p>
i) combined.loom - Loom hd5 file containing the combined data for all the developmental timepoints used in the publication.
</p>
ii) wt.loom - Loom hd5 file containing the processed data for WT mice.
</p>
iii) mafko.loom - Loom hd5 file containing the processed data for Maf KO mice.
</p>
iv) development.rds - .RDS file containing the processed Seurat object for all developmental timepoints
</p>
v) mafko.rds - .RDS file containing the processed Seurat object for combineed WT and Maf KO mice scRNA-seq data.

### Contact
Any questions, feel free to post on the Issues section. If you require any of the processed files reach out to:
jmg2008@med.cornell.edu or fri2002@med.cornell.edu
