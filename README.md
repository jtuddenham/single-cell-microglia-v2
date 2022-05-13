# single-cell-microglia-v2
Code for the analysis and figure generation for the manuscript "A cross-disease human microglial framework identifies disease-enriched subsets and tool compounds for microglial polarization."

raw_data: contains h5 output files for all datasets included in this study. Mapping between sequencing ID and donor ID can be found in Table S1 in intermediate_data.
intermediate_data: contains example outputs from our clustering/diffex pipeline, such as differentially expressed genes for all clusters and final cluster labels.
Preprocessing and Clustering: contains scripts and helper functions for clustering/differential expression, validation (bootstrap and constellation/reclassification), and label transfer. 

Dataset used for label transfer is from:
Hasselmann J, et al. Development of a Chimeric Model to Study and Manipulate Human Microglia In Vivo. Neuron. 2019 Sep 25;103(6):1016-1033.e10. doi: 10.1016/j.neuron.2019.07.002. Epub 2019 Jul 30. PMID: 31375314; PMCID: PMC7138101.
Available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133432
