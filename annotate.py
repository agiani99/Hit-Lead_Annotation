"""Annotate a ChEMBL activity CSV and build the radar-plot HTML report.

    python annotate.py CHEMBL2808_EC50.csv
"""
import argparse
from pathlib import Path

from hitanno import pipeline, report, surechembl

here = Path(__file__).parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
ap.add_argument("input", help="CSV with SMILES and pChEMBL columns (',' ';' or tab separated)")
ap.add_argument("--out-dir", default=str(here / "output"))
ap.add_argument("--config", default=str(here / "config.yaml"), help="axes + thresholds")
ap.add_argument("--models-dir", help="folder with your own ADMET-AI-format models (default: bundled ADMET-AI)")
ap.add_argument("--surechembl", default=surechembl.REMOTE,
                help="SureChEMBL bulk_data folder: URL (default, streamed) or local directory with the parquet files")
ap.add_argument("--patent-classes", default=surechembl.DEFAULT_CLASSES,
                help="regex an IPC/CPC class must match for a patent to count ('' = any patent); filters text-mining noise")
ap.add_argument("--patents-min-phase", type=int, default=0, help="only look up patents for compounds at/above this phase (faster)")
ap.add_argument("--patent-term", type=int, default=25, help="years after the first patent date after which it is shown as expired (green)")
ap.add_argument("--no-surechembl", action="store_true")
ap.add_argument("--no-admet", action="store_true")
ap.add_argument("--min-phase", type=int, default=1, help="lowest max phase plotted (0 is never plotted)")
ap.add_argument("--batch-size", type=int, default=40, help="ChEMBL molecules per request")
ap.add_argument("--cache-dir", default=str(here / "cache"))
ap.add_argument("--limit", type=int, help="only read the first N input rows (testing)")
ap.add_argument("--plots-only", action="store_true", help="rebuild HTML from an existing output CSV")
a = ap.parse_args()

out = Path(a.out_dir)
out.mkdir(parents=True, exist_ok=True)
cfg = pipeline.load_config(a.config)
stem = Path(a.input).stem
csv_out = out / f"{stem}_annotated.csv"

if a.plots_only:
    import pandas as pd
    df = pd.read_csv(csv_out)
else:
    raw = pipeline.read_input(a.input, a.limit)
    df = pipeline.annotate(cfg, raw, a.cache_dir, a.models_dir, a.surechembl, not a.no_surechembl,
                           not a.no_admet, a.batch_size, a.patent_classes, a.patents_min_phase)
    df.to_csv(csv_out, index=False)
    print(f"Wrote {csv_out}")
report.build_report(df, cfg, out / f"{stem}_radar.html", max(a.min_phase, 1), patent_term=a.patent_term)
print(f"Wrote {out / (stem + '_radar.html')}")
