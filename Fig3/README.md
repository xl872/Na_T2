# Fig. 3 and Supplementary Tables 1-2

| Panel | Script | Input data | Output |
|---|---|---|---|
| 3B | none | `data/DWI_scan613`, `data/DWI_scan627` | none |
| 3C, ADC maps | `scripts/myADCProt_scan613.m`, `scripts/myADCProt_scan627.m` | `data/DWI_scan613`, `data/DWI_scan627`, `data/bruker` | `ADC+orig`, `ADC_R2+orig` |
| 3C, Supplementary Table 1 | `scripts/ADCstrip_Fig3C.m` | `data/DWI_scan613`, `data/DWI_scan627`, `data/bruker` | `results/ADCstrip_Fig3C_cells.csv`, `results/ADCstrip_Fig3C_decay.csv`, figure |
| 3D, Supplementary Table 2 | `scripts/Fig3D_ADCfastslow.m` | `data/Fig3D/Fig3D_input_signals.csv`, `data/Fig3D/Fig3D_input_sheetADC.csv` | `results/Fig3D_ADC_table.csv`, figures |

To run, from the repository root: `cd Fig3/data && unzip -oq bruker_parameters.zip && for d in DWI_*; do (cd "$d" && unzip -oq "${d}_afni.zip"); done`, then in MATLAB run each script from `Fig3/scripts`, top row first.
