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
