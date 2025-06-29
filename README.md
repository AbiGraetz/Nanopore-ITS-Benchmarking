# Nanopore-ITS-Benchmarking
This repo contains scripts and intermediate files used for benchmarking analysis used in:
Graetz and Feng et al., 2025, 'Benchmarking fungal species classification using Oxford Nanopore Technologies long-read metabarcodes'

There are eight directories: 
1. qc-pipeline: Scripts and input files for quality control of fungal mock community sequences. Written by A. Ringeri.
2. JF_code: Scripts and input/intermediate files for sequence classification, generating F1 score dotplot figures, and machine learning overfitting analyses. Written by J. Feng.
3. SimulateNanoporeMetabarcodes: Script used to create the simNCBI database for training machine learning classifiers. Written by A. Graetz.
4. RandomForestRegression: Script and input files for random forest regression simulations of F1 score variance. Written by A. Graetz. 
5. MinimumSequencingDepth: Script and input files for simulating minimum required sequences to detect all mock community species at minimum abundance thresholds. Written by A. Graetz. 
6. UnseenSpecies: Script and input files for assessing classifier behaviour with incomplete reference databases. Written by A. Graetz. 
7. Mock3_CandidaDiagnostics: Script and input files for diagnostic simulation using *Candida tropicalis*. Written by A. Graetz. 
8. AlphaDiversity_BhattacharyyaDistance: Script and input files for community diversity analyses. Written by A. Graetz. 

While input files are provided, raw data can be accessed from Zenodo: 10.5281/zenodo.15751794