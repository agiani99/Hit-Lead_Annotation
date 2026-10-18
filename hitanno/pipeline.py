"""Annotation pipeline: CSV in -> one annotated row per unique compound."""
import csv
import math
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from . import admet, surechembl
from .chembl import Chembl

RDLogger.DisableLog("rdApp.*")

COLS = {  # accepted (lower-case) header names
    "smiles": ["smiles", "canonical_smiles"],
    "pchembl": ["pchembl value", "pchembl_value", "pchembl"],
    "id": ["molecule chembl id", "molecule_chembl_id", "chembl_id", "chembl id"],
    "name": ["molecule name", "name", "pref_name"],
    "phase": ["molecule max phase", "max_phase"],
}
RDKIT_PROPS = {
    "alogp": Crippen.MolLogP, "psa": rdMolDescriptors.CalcTPSA, "rtb": rdMolDescriptors.CalcNumRotatableBonds,
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(), "hbd": rdMolDescriptors.CalcNumHBD,
}


def load_config(path):
    return yaml.safe_load(Path(path).read_text())["axes"]


def read_input(path, limit=None):
    with open(path, newline="", encoding="utf-8-sig") as f:
        sep = csv.Sniffer().sniff(f.read(4096), delimiters=";,\t").delimiter
    df = pd.read_csv(path, sep=sep, dtype=str)
    low = {c.lower().strip(): c for c in df.columns}
    pick = {k: next((low[n] for n in names if n in low), None) for k, names in COLS.items()}
    if not pick["smiles"]:
        raise SystemExit(f"No SMILES column found in {path} (columns: {list(df.columns)})")
    out = pd.DataFrame({k: df[v] if v else None for k, v in pick.items()})
    out["pchembl"] = pd.to_numeric(out["pchembl"], errors="coerce")
    out["phase"] = pd.to_numeric(out["phase"], errors="coerce")
    return out.head(limit) if limit else out


_uncharger = rdMolStandardize.Uncharger()
_chooser = rdMolStandardize.LargestFragmentChooser()


def standardize(smiles):
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None, None, None
    m = _uncharger.uncharge(_chooser.choose(m))
    return Chem.MolToSmiles(m), Chem.MolToInchiKey(m), m


def unique_compounds(raw, log=print):
    rows = []
    for r in raw.itertuples(index=False):
        smi, ik, mol = standardize(r.smiles)
        if ik is None:
            log(f"  skipped unparsable SMILES: {r.smiles!r}")
            continue
        rows.append({"inchikey": ik, "smiles": smi, "chembl_id": r.id, "name": r.name, "csv_phase": r.phase,
                     "pchembl_input": r.pchembl})
    df = pd.DataFrame(rows)
    g = df.groupby("inchikey", sort=False)
    out = g.agg(smiles=("smiles", "first"),
                chembl_id=("chembl_id", lambda s: s.dropna().iloc[0] if s.notna().any() else None),
                name=("name", lambda s: s.dropna().iloc[0] if s.notna().any() else None),
                csv_phase=("csv_phase", "max"),
                pchembl=("pchembl_input", "max"),  # highest potency across assays
                n_input_rows=("smiles", "size")).reset_index()
    return out


def _strict(values, direction):
    """Strictest measurement: highest for liabilities, lowest when higher is better."""
    return min(values) if direction == "higher_better" else max(values)


def annotate(cfg, raw, cache_dir, models_dir=None, surechembl_src=surechembl.REMOTE, use_surechembl=True,
             use_admet=True, batch_size=40, patent_classes=surechembl.DEFAULT_CLASSES,
             patents_min_phase=0, log=print):
    log(f"Input rows: {len(raw)}")
    df = unique_compounds(raw, log)
    log(f"Unique compounds (by InChIKey): {len(df)}")
    ch = Chembl(cache_dir, batch_size)

    # --- resolve ChEMBL IDs for rows lacking one
    for i in df.index[df.chembl_id.isna()]:
        df.at[i, "chembl_id"] = ch.id_from_inchikey(df.at[i, "inchikey"])
    ids = df.chembl_id.dropna().tolist()

    log("ChEMBL: molecules / phase / properties")
    mols = ch.molecules(ids)
    df["name"] = [mols.get(c, {}).get("name") or n or c for c, n in zip(df.chembl_id, df["name"])]
    df["max_phase"] = [mols[c]["max_phase"] if c in mols else (0 if math.isnan(p) else int(p))
                       for c, p in zip(df.chembl_id, df.csv_phase)]

    log("ChEMBL: measured lighthouse activities")
    act = {}
    for ax in cfg:
        if ax["kind"] == "activity":
            act[ax["key"]] = ch.pchembl_by_target(ids, ax["targets"])
        elif ax["kind"] == "cytotox":
            act[ax["key"]] = ch.pchembl_by_target(ids, None, standard_type="CC50")
    mw = {c: float(mols[c]["props"].get("full_mwt") or 0) for c in ids if c in mols}
    logs = ch.logs_values(mw) if any(a["kind"] == "solub" for a in cfg) else {}

    pred = pd.DataFrame()
    if use_admet:
        log("ADMET-AI predictions")
        pred = admet.predict(df.smiles.tolist(), models_dir, log=log)

    rdmols = {s: Chem.MolFromSmiles(s) for s in df.smiles}
    for ax in cfg:
        k, vals, srcs = ax["key"], [], []
        for r in df.itertuples():
            v, s = None, ""
            if ax["kind"] == "prop":
                p = mols.get(r.chembl_id, {}).get("props", {})
                if p.get(ax["field"]) not in (None, ""):
                    v, s = float(p[ax["field"]]), "chembl"
                else:
                    v, s = RDKIT_PROPS[ax["field"]](rdmols[r.smiles]), "rdkit"
            else:
                measured = []
                if ax["kind"] in ("activity", "cytotox"):
                    measured = [x for lst in act[k].get(r.chembl_id, {}).values() for x in lst]
                elif ax["kind"] == "solub":
                    measured = logs.get(r.chembl_id, [])
                if measured:
                    v, s = _strict(measured, ax["direction"]), "measured"
                elif ax.get("pred") in pred.columns:
                    pv = pred.at[r.smiles, ax["pred"]] if r.smiles in pred.index else None
                    if pv is not None and not pd.isna(pv):
                        v, s = float(pv), "predicted"
            vals.append(v)
            srcs.append(s)
        df[k], df[f"{k}_source"] = vals, srcs

    if use_surechembl:
        log("SureChEMBL: first patent per compound")
        todo = df[df.max_phase >= patents_min_phase]
        pat = surechembl.first_patents(todo.inchikey.tolist(), surechembl_src, class_regex=patent_classes, log=log)
        df = df.merge(pat, on="inchikey", how="left")
    else:
        for c in ["patent_number", "patent_pub_date", "patent_originator", "patent_use", "patent_ipc", "patent_match"]:
            df[c] = None
    return df.drop(columns=["csv_phase"])
