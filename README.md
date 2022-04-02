# "Transcriptional specification of fetal to adult liver  vascular functional zonation"
Code for reproducibility and downstream analysis of the single cell RNA sequencing data available in "Transcriptional specification of fetal to adult liver  vascular functional zonation" [[add link/doi]], including: i) Single cell RNA sequencing for developmental timepoints in mouse and ii) Single cell RNA sequencing of WT and Maf KO endothelial cells

## Scripts
This repository contains the scripts used for analysis of single cell RNA seq data for "Transcriptional specification of fetal to adult liver  vascular functional zonation" [[add link/doi]], including the following files:
</p>
i) <strong>development.R</strong> - Seurat base analysis of single cell RNA seq data for developmental timepoints.
</p>
ii) <strong>scVelo_development.ipynb</strong>  RNA velocity analysis using velocyto (https://github.com/velocyto-team/velocyto.py) and scVelo (https://github.com/theislab/scvelo) for developmental timepoints to infer transitions and genes driving the transcriptional changes. The .ipynb file contains the Jupyter notebook and includes graphical output.
</p>
iii) <strong>MafKO.R</strong> Seurat base analysis for single cell RNA seq data for WT vs Maf KO experiment
</p>
iv) <strong>scVelo_MafKO.ipynb</strong> RNA velocity analysis for either WT and Maf knockout mice, for comparison of developmental trajectories of endothelial cells upon Maf knockout. The .ipynb file contains the Jupyter notebook and includes graphical output.

## Metadata files
We included the metadata files containing cell barcode, QC data, clustering, cluster names and UMAP coordinates for both development and MafKO datasets as .csv files:
</p>
i) <strong>development_metadata.csv</strong>
</p>
ii) <strong>MafKO_metadata.csv</strong>
</p>
These files are also available through the GEO repository (GSE174209 and GSE174208).

## Processed files available at GEO
Data for this publication is available at GEO (GSE174209 and GSE174208). In addition to raw data, we provide the Seurat objects for both development and c-Maf KO  to facilitate downstream analysis without the need to re-process the raw data. 

### Contact
Any questions, feel free to post on the Issues section. If you require any of the processed files reach out to:
jmg2008@med.cornell.edu or fri2002@med.cornell.edu
