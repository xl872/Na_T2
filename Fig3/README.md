# Fig. 3 and Supplementary Tables 3-5

| Panel | Script | Input data | Output |
|---|---|---|---|
| 3A | `scripts/T2estimate_M1T2.m`, `scripts/T2estimate_M3T2.m` | `data/M1`, `data/M3` (repository root) | `T2`, `T2short`, `T2long`, `R2bi` maps |
| Supplementary Tables 4-5, #3 and #4 | `scripts/T2estimate_metrics_M1T2.m`, `scripts/T2estimate_metrics_M3T2.m` | `data/M1`, `data/M3` | maps and fits for those tables |
| Supplementary Tables 3-5, #3 and #4 | `scripts/R2mono_table.m`, `scripts/T2map_ROIcontour_p2.m` | `data/M1`, `data/M3` | `results/T2_ROI_stats_p2.csv`, map figures |
| 3B | `scripts/SNRfit_Fig3B.m`, `scripts/R2bi_ROIbar_fig3.m` | `data/M1`, `data/M3` | plots |
| Supplementary Tables 4-5, #5 | `scripts/T2estimate_metrics_scan142_160.m` | `data/Na23_scan142-160` | `results/T2_ROI_stats_scan142-160.csv` |
| Supplementary Tables 4-5, #6 | `scripts/T2estimate_metrics_scan273_293.m` | `data/Na23_scan273-293`, `data/bruker` | `results/T2_ROI_stats_scan273-293.csv` |
| 3C | none | `data/Fig3C_values.csv` | plotted values |
| 3A, scripts at the repository root | `scripts/T2estimate_M1.m`, `scripts/T2estimate_M3.m` | `data/M1`, `data/M3` | `T2`, `T2short`, `T2long` maps |

To run, from the repository root: `(cd data/M1 && unzip -oq 2025051802Na.M1T2.rawdata.zip && unzip -oq ../../Fig3/data/M1T2_additional.zip); (cd data/M3 && unzip -oq 20250518Na.M3T2.rawdata.zip && unzip -oq ../../Fig3/data/M3T2_additional.zip); cd Fig3/data && unzip -oq bruker_parameters.zip && for d in Na23_*; do (cd "$d" && unzip -oq "${d}_afni.zip"); done`, then in MATLAB run each script from `Fig3/scripts`, top row first.
