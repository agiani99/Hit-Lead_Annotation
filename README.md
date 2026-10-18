# HIT_anno

Annotates a ChEMBL activity CSV (SMILES + pChEMBL) with patent, property, ADMET and development-phase data,
and builds a static HTML page of radar plots (phase 4 -> 1, three per row).

```
pip install -r requirements.txt
python annotate.py CHEMBL2808_EC50.csv --patents-min-phase 1
```

Output in `output/`: `<name>_annotated.csv` and `<name>_radar.html`. ChEMBL calls are cached in `cache/`.

## Options
| flag | meaning |
|---|---|
| `--config FILE` | axes and cutoffs (`config.yaml`: `thr` measured, `pred_thr` predicted) |
| `--plots-only` | re-render the HTML from the existing CSV, e.g. after editing cutoffs |
| `--models-dir DIR` | your own ADMET-AI-format models instead of the bundled ones |
| `--patents-min-phase N` | patent lookup only for phase >= N (much faster) |
| `--patent-classes REGEX` | IPC/CPC filter for the first patent (default `^(C07\|A61K\|A61P)`, `''` = any) |
| `--patent-term 25` | years after which the patent date is shown green (expired) |
| `--surechembl DIR` | local SureChEMBL parquet files (default: streamed from EBI) |
| `--min-phase N` | lowest phase plotted (minimum 1) |

## Rules
* One row per InChIKey; potency = highest pChEMBL.
* One value per axis: measured (ChEMBL) beats predicted (ADMET-AI); the strictest measurement is kept.
* Radar: centre good, green disc = allowed, red `!` = violation, black dot = measured, blue dot = predicted.

## Limits
* Patent date is the publication date, not the application date; "first patent" comes from noisy text mining, so check it.
* No AUC or cytotoxicity predictor by default (axes appear only with measured data or a local model).
* PubChem is not queried.
