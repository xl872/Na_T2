# Fig. 5 and Supplementary Table 6

| Panel | Script | Input data | Output |
|---|---|---|---|
| Maps for 5A-D | `scripts/T2estimate_scan142_160.m`, `scripts/T2estimate_scan402_418.m`, `scripts/myADCProt_scan175.m`, `scripts/myADCProt_scan613.m` | `data/Na23_scan142-160`, `data/Na23_scan402-418`, `data/DWI_scan175`, `data/DWI_scan613`, `data/bruker` | T2* and ADC maps |
| 5C-D, Supplementary Table 6 | `scripts/alignT2toDTI175.m`, `scripts/ROI_T2long_T2vsADC_scan175.m` | `data/DWI_scan175`, `data/Na23_scan142-160`, `data/bruker` | `results/ROI_T2long_*_scan175.csv`, figures |
| 5A-B, Supplementary Table 6 | `scripts/alignT2toDTI613New.m`, `scripts/ROI_T2long_T2vsADC_scan613.m` | `data/DWI_scan613`, `data/Na23_scan402-418`, `data/bruker` | `results/ROI_T2long_*_scan613.csv`, figures |

To run, from the repository root: `cd Fig5/data && unzip -oq bruker_parameters.zip && for d in DWI_* Na23_*; do (cd "$d" && unzip -oq "${d}_afni.zip"); done`, then in MATLAB run each script from `Fig5/scripts`, top row first.
