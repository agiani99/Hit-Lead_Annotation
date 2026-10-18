"""ChEMBL REST access with on-disk caching and batching."""
import hashlib
import json
import time
from pathlib import Path

import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"
UNIT_TO_G_PER_L = {"ug.mL-1": 1e-3, "mg.mL-1": 1.0, "ug/ml": 1e-3, "mg/ml": 1.0, "g.L-1": 1.0, "mg.L-1": 1e-3}
UNIT_TO_MOL_PER_L = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9}


class Chembl:
    def __init__(self, cache_dir, batch_size=40):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.batch = batch_size
        self.http = requests.Session()

    def _get(self, path, params):
        key = hashlib.md5((path + json.dumps(params, sort_keys=True)).encode()).hexdigest()
        f = self.cache / f"{key}.json"
        if f.exists():
            return json.loads(f.read_text())
        for attempt in range(5):
            try:
                r = self.http.get(f"{BASE}/{path}.json", params=params, timeout=90)
                if r.status_code == 200:
                    f.write_text(r.text)
                    return r.json()
                if r.status_code in (404, 400):
                    return {}
            except requests.RequestException:
                pass
            time.sleep(2 ** attempt)
        raise RuntimeError(f"ChEMBL request failed: {path} {params}")

    def _paged(self, path, params, key):
        out, offset = [], 0
        while True:
            d = self._get(path, {**params, "limit": 1000, "offset": offset})
            out += d.get(key, [])
            if not d or not d.get("page_meta", {}).get("next"):
                return out
            offset += 1000

    def molecules(self, chembl_ids):
        """chembl_id -> {name, max_phase, props{}, inchikey}"""
        ids, res = sorted(set(chembl_ids)), {}
        for i in range(0, len(ids), self.batch):
            chunk = ids[i:i + self.batch]
            for m in self._paged("molecule", {"molecule_chembl_id__in": ",".join(chunk)}, "molecules"):
                try:
                    phase = int(float(m.get("max_phase")))
                except (TypeError, ValueError):
                    phase = 0
                res[m["molecule_chembl_id"]] = {
                    "name": m.get("pref_name"),
                    "max_phase": phase,
                    "props": m.get("molecule_properties") or {},
                    "inchikey": (m.get("molecule_structures") or {}).get("standard_inchi_key"),
                }
        return res

    def id_from_inchikey(self, inchikey):
        ms = self._get("molecule", {"molecule_structures__standard_inchi_key": inchikey}).get("molecules", [])
        return ms[0]["molecule_chembl_id"] if ms else None

    def pchembl_by_target(self, chembl_ids, targets, standard_type=None):
        """{chembl_id: {target: [pchembl,...]}} for pchembl-bearing activities."""
        ids, res = sorted(set(chembl_ids)), {}
        for i in range(0, len(ids), self.batch):
            p = {"molecule_chembl_id__in": ",".join(ids[i:i + self.batch]), "pchembl_value__isnull": "false"}
            if targets:
                p["target_chembl_id__in"] = ",".join(targets)
            if standard_type:
                p["standard_type"] = standard_type
            for a in self._paged("activity", p, "activities"):
                res.setdefault(a["molecule_chembl_id"], {}).setdefault(a["target_chembl_id"], []).append(
                    float(a["pchembl_value"]))
        return res

    def logs_values(self, mols):
        """mols: {chembl_id: molecular_weight}. Returns {chembl_id: [logS,...]} from Solubility activities."""
        import math
        ids, res = sorted(mols), {}
        for i in range(0, len(ids), self.batch):
            p = {"molecule_chembl_id__in": ",".join(ids[i:i + self.batch]), "standard_type": "Solubility"}
            for a in self._paged("activity", p, "activities"):
                try:
                    v, u, mw = float(a["standard_value"]), a["standard_units"], mols[a["molecule_chembl_id"]]
                except (TypeError, ValueError, KeyError):
                    continue
                if v <= 0 or a.get("standard_relation") not in ("=", None):
                    continue
                if u in UNIT_TO_MOL_PER_L:
                    mol_l = v * UNIT_TO_MOL_PER_L[u]
                elif u in UNIT_TO_G_PER_L and mw:
                    mol_l = v * UNIT_TO_G_PER_L[u] / mw
                else:
                    continue
                res.setdefault(a["molecule_chembl_id"], []).append(math.log10(mol_l))
        return res
