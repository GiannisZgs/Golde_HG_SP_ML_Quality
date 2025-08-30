# Systematic Benchmarking of a Noise-tolerant Conductive Hydrogel electrode for Epidermal Bioelectronics

This repository contains the codebase for analyzing and comparing ECG signals recorded with conventional AgCl electrodes and polymer-printed hydrogel (PPHG) electrodes. The analysis pipeline processes raw ECG data, extracts features, performs clustering analysis, and generates visualizations for quality comparison. Additionally, the comparison is extended to other biosignals modalities, EEG, EMG and EOG.

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

### EEG Analysis

Run the following scripts to analyze EEG data with and without preprocessing:

```matlab
cd scripts/other_modalities
```

1. Analysis with preprocessing:
```matlab
run('EEG_with_preprocessing.m')
```
Performs full EEG analysis with bandpass filter (0.5–40 Hz) and detrending applied directly in MATLAB code. The workflow includes PSD computation, spectrogram and CWT visualization, segmentation into baseline / thinking / eyes-closed conditions, delta power extraction (0.5–4 Hz), statistical comparison using t-tests, and time-resolved power plots with shaded intervals for each condition. The results provided from this script were imported into OriginLab for plotting.
### EMG Analysis

2. Analysis without preprocessing:
```matlab
run('EEG_without_preprocessing.m')
```
Implements the same analysis pipeline as EEG_with_preprocessing.m, except that filtering and detrending are performed beforehand using MATLAB’s Signal Analyzer GUI. Used to compare results from GUI-preprocessed data versus in-code preprocessing. The results provided from this script were imported into OriginLab for plotting.


### EMG Analysis

```matlab
run('EMG_with_preprocessing.m')
```

Reads EMG data from AgCl and hydrogel electrodes, applies preprocessing via the helper function preprocess_emg_complete_adv.m, and then performs envelope extraction (Hilbert transform), spectrogram and CWT plotting, and peak-to-peak voltage analysis with comparative bar charts. The results provided from this script were imported into OriginLab for plotting.

### EOG Analysis

```matlab
run('EOG_with_preprocessing_and_analysis_complete_code 1.m')
```

Processes raw EOG recordings by applying multiple cleaning steps (filtering, notch, baseline correction, detrending, wavelet denoising, median filtering). The cleaned data are segmented into defined activities (“look up”, “blink 3x”, “look down”), from which maximum, minimum, and peak-to-peak values are extracted. Statistical analysis includes one-way ANOVA with Tukey post-hoc tests, and visual outputs include preprocessing stage plots and bar graphs for amplitude comparisons. The results provided from this script were imported into OriginLab for plotting.



