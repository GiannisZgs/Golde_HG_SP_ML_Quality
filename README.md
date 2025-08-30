# Systematic Benchmarking of a Noise-tolerant Conductive Hydrogel electrode for Epidermal Bioelectronics

This repository contains the codebase for analyzing and comparing ECG signals recorded with conventional AgCl electrodes and polymer-printed hydrogel (PPHG) electrodes. The analysis pipeline processes raw ECG data, extracts features, performs clustering analysis, and generates visualizations for quality comparison.

## Repository Structure

```
Golde_HG_SP_ML_Quality/
├── data/                  # Output data from MATLAB processing
├── imgs_figures/          # Generated visualizations
├── scripts/               # MATLAB processing scripts
│   ├── clustering/        # Clustering analysis scripts
│   ├── ECG_processing/    # ECG processing scripts
│   └── other_modalities/  # Scripts for EEG, EMG, EOG analysis
├── utils/                 # Utility functions
│   ├── clustering/        # Clustering utility functions
│   ├── ECG_preprocessing/ # ECG preprocessing functions
│   ├── feature_extraction/# Feature extraction functions
│   └── other/             # Miscellaneous utilities
└── visualization/         # R scripts for visualization
```

## Prerequisites

- MATLAB R2020b or newer with Signal Processing Toolbox
- Python 3.7+ with scipy, numpy
- R 4.0+ with required packages:
  - jsonlite, ggplot2, dplyr, tidyr, patchwork, viridis, scales, grid, cowplot, ggrepel, ggbeeswarm

## Workflow to Reproduce Results

### 1. Set Up Environment

1. Clone this repository
2. Run MATLAB and navigate to the repository root folder
3. Run the setup script to add all necessary directories to the MATLAB path:
   ```matlab
   run('utils/setup_environment.m')
   ```

   Alternatively, this will be automatically be handled by executing any of the below ECG processing scripts (see ##2.).

### 2. Process ECG Data (MATLAB)

The processing pipeline consists of several steps. Please check each file and adjust 
the boolean flags at the start of each script to control plotting, saving of files and overall behavior:

1. Extract and preprocess ECG data:
   ```matlab
   cd scripts/ECG_processing
   run('extract_all_ecgs.m')
   ```
   *: For the manually cleaned recordings experiment, run this instead:
    run('extract_manually_cleaned_ecgs_p1_p5_p10_p39.m')

2. Perform ECG analysis to extract heartbeat profiles:
   ```matlab
   run('ECG_analysis.m')
   ```

3. Run clustering analysis for participant identification:
   ```matlab
   cd ../clustering
   run('clustering_analysis_participants.m')
   ```

4. Run clustering analysis for channel identification:
   ```matlab
   run('clustering_analysis_channels.m')
   ```

These scripts will generate `.mat` files in the data directory.

### 3. Convert MATLAB Outputs to JSON

Use the Python converter to transform MATLAB outputs to JSON format:

```bash
cd visualization
python convert_mat_to_json.py
```

This will convert the `.mat` files to `.json` files in the data directory.

### 4. Generate Visualizations (R)

Run the R scripts to generate visualizations from the JSON data:

```bash
cd visualization
Rscript figure_5a_S9.R     # ECG waveform overlay visualization
Rscript figure_5b_S10.R    # Channel clustering visualization
Rscript figure_5c_S11.R    # Participant clustering visualization
Rscript figure_5d_S12.R    # ECG similarity metric visualization
Rscript figure_5e.R        # ECG identification performance visualization
Rscript figure_4c_S6.R     # ECG signal processing visualization
Rscript figure_4d.R        # ECG feature annotation visualization
Rscript figure_4e_S8.R     # ECG noise reduction visualization
Rscript figure_S7.R        # Power spectral density visualization
```

Visualizations will be saved in the `imgs_figures` directory.

## Output Directories

- `data/`: Contains the processed data outputs from MATLAB (.mat and .json files)  
- `imgs_figures/`: Contains visualization outputs from R scripts

## Key Results

The main outcomes of the ECG analysis pipeline:
- Comparative quality assessment between AgCl and PPHG electrodes
- Participant identification performance 
- Channel (lead) identification performance
- Noise reduction and signal quality metrics across electrode types
- Feature extraction and comparison between electrode types

## Additional Modalities Analysis

Besides the ECG analysis pipeline, this repository also includes scripts for analyzing other bio-signal modalities:

### EEG Analysis

Run the following scripts to analyze EEG data with and without preprocessing:

```matlab
cd scripts/other_modalities
```
1. Analysis without preprocessing:
```matlab
run('EEG_without_preprocessing.m')
```
This script performs frequency domain analysis of raw EEG data without filtering or artifact removal, focusing on delta band power differences between eyes-open and eyes-closed conditions.

2. Analysis with preprocessing:
```matlab
run('EEG_with_preprocessing.m')
```
This script applies bandpass filtering and artifact removal to EEG data before analysis, allowing comparison of the impact of preprocessing on spectral power differences.

### EMG Analysis

For EMG signal processing and visualization:

```matlab
run('EMG_with_preprocessing.m')
```

This script loads, preprocesses, and compares EMG signals from AgCl and hydrogel electrodes, generating:
- Time-domain visualizations with signal envelopes
- Frequency-domain analysis
- Time-frequency representations (spectrograms and scalograms)
- Peak-to-peak voltage comparisons between electrode types




