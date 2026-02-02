# FLKeys-Coral-Seascape-Gen
Scipts and data to accompany paper:\
​### Black, Kristina L., John P. Rippe, and Mikhail V. Matz. 2025. “Environmental Drivers of Genetic Adaptation in Two Corals from the Florida Keys.” Evolutionary Applications, 18:e70126

## Key datasets:
- **Agaricia_nozoox.ibsMat**: matrix of genetic distances for *A. agaricites* (identity-by-state), computed from *de novo* 2b-RAD data by ANGSD.
- **Agaricia_env_feb2025.RData**: Metadata for *A. agaricites*: environmental data (interpolateed from SERC and ERDDAP) and geographical coordinates of each sample.
- **Agaricia_admix.csv**: by-lineage admixture proportions
- **Porites_nozoox.ibsMat**: matrix of genetic distances for *P.astreoides* (identity-by-state), computed from *de novo* 2b-RAD data by ANGSD.
- **Porites_env_feb2025.RData**: Metadata for *P.astreoides*: environmental data (interpolateed from SERC and ERDDAP) and lon, lat of each sample.
- **Porites_admix.csv**: by-lineage admixture proportions
- **SERC_krigs_bottom_prepost_5k.RData**: environmental data for 5km raster of Florida Keys Reef Tract, for 2003-2007 ("pre") and 2015-2019 ("post")

## Key scripts:
- **Agaricia_dec2024_yellow.R** (also, similarly scripts for *Porites* and other lineages including all lineages together, "all"): RDAforest analysis and predictive model building for a specific dataset.

