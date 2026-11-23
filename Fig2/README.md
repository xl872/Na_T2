# Fig. 2 and Supplementary Tables 1-2

| Panel | Script | Input data | Output |
|---|---|---|---|
| 2B | none | `data/DWI_scan613`, `data/DWI_scan627` | none |
| 2C, ADC maps | `scripts/myADCProt_scan613.m`, `scripts/myADCProt_scan627.m` | `data/DWI_scan613`, `data/DWI_scan627`, `data/bruker` | `ADC+orig`, `ADC_R2+orig` |
| 2C, Supplementary Table 1 | `scripts/ADCstrip_Fig2C.m` | `data/DWI_scan613`, `data/DWI_scan627`, `data/bruker` | `results/ADCstrip_Fig2C_cells.csv`, `results/ADCstrip_Fig2C_decay.csv`, figure |
| 2D, Supplementary Table 2 | `scripts/Fig2D_ADCfastslow.m` | `data/Fig2D/Fig2D_input_signals.csv`, `data/Fig2D/Fig2D_input_sheetADC.csv` | `results/Fig2D_ADC_table.csv`, figures |

To run, from the repository root: `cd Fig2/data && unzip -oq bruker_parameters.zip && for d in DWI_*; do (cd "$d" && unzip -oq "${d}_afni.zip"); done`, then in MATLAB run each script from `Fig2/scripts`, top row first.
