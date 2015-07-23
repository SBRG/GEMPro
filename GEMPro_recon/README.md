# GEMPro_recon

In this directory you will find the main dataframes and structure files for a given model. The directory structure is as follows:

$MODEL
├── data_frames
├── model_files
├── structure_files
│   ├── experimental  # all available pdb structures matched to the model (unfiltered!)
│   ├── homology_models  # all available homology models matched to the model (chain added & cleaned where possible)
│   ├── ssb_best  # all files matching the ssb_best_file column in the main dataframe
└── tutorials  # any ipython tutorials for using GEM-PRO
