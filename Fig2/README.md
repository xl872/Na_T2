# Fig. 2

| Panel | Script | Input data | Output |
|---|---|---|---|
| 2A | none | `data/Na23_scan142-160`, `data/DWI_scan175` | none |
| 2B, ADC map | `scripts/myADCProt_scan175.m` | `data/DWI_scan175`, `data/bruker` | `ADC+orig`, `ADC_R2+orig` |
| 2B, curve | none | `data/Fig2B_ROI_signal.csv` | source data for the plotted curve |
| 2C | `scripts/T2estimate_scan142_160.m` | `data/Na23_scan142-160` | `T2f`, `T2s`, `T2` maps, decay fit figure |

To run, from the repository root: `cd Fig2/data && unzip -oq bruker_parameters.zip && for d in DWI_* Na23_*; do (cd "$d" && unzip -oq "${d}_afni.zip"); done`, then in MATLAB run each script from `Fig2/scripts`, top row first.
