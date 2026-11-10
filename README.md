Yu G, Liu X, ... Zhou XA. *High-Field Multinuclear MRI Reveals Sodium Relaxation Heterogeneity in Fixed Cortical Organoids.* Manuscript submitted.

## Folders

- `Fig2/`: scripts and data for Fig. 2
- `Fig3/`: scripts, data and results for Fig. 3 and Supplementary Tables 1-2
- `Fig4/`: scripts, data and results for Fig. 4 and Supplementary Tables 3-5
- `Fig5/`: scripts, data and results for Fig. 5 and Supplementary Table 6
- `conversion/`: conversion of the Bruker data to AFNI format
- `helper/`: shared MATLAB function

The existing `scripts/T2estimate_M1.m`, `scripts/T2estimate_M3.m`, `data/M1` and `data/M3` stay where they are; Fig. 4 reads `data/M1` and `data/M3`.

## Software

- MATLAB R2025a with the Curve Fitting, Image Processing and Statistics toolboxes
- AFNI 24.0.09
- afni_matlab, for BrikLoad and WriteBrik

## Data availability

Raw data: Harvard Dataverse, https://doi.org/10.7910/DVN/HRVBRW

# T2 Estimation

MATLAB script for estimating T2 from multi-TE EPI data using mono-exponential and bi-exponential fitting.

## Main script
- Loads multi-TE AFNI BRIK datasets
- Estimates background noise from edge voxels
- Applies light 3D Gaussian smoothing
- Computes ROI-averaged signal decay across TE
- Fits:
  - single-exponential model for T2
  - bi-exponential model for short/long T2 components
- Performs voxel-wise fitting to generate T2 maps
- Saves AFNI output volumes and summary figures

## Required inputs

- Multi-TE AFNI datasets  
  `2025051802Na.<run>.an200epi+orig.BRIK`
- `ROI+orig.BRIK`
- `formask+orig.BRIK`

## Outputs

- ROI fitting results (`.mat`)
- T2 maps:
  - `T2`
  - `T2short`
  - `T2long`
- Masked T2 maps
- QC and fitting figures (`.jpg`, `.eps`)

## Requirements

- MATLAB
- Image Processing Toolbox
- Curve Fitting Toolbox
- AFNI MATLAB functions:
  - `BrikLoad`
  - `WriteBrik`
