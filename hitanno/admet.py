"""ADMET predictions via ADMET-AI (optionally with a local models folder)."""
import contextlib
import io
import os

import pandas as pd


def predict(smiles, models_dir=None, chunk=500, log=print):
    """Return DataFrame indexed by SMILES with all model outputs. Batched to bound memory."""
    from admet_ai import ADMETModel
    kw = {"models_dir": models_dir} if models_dir else {}
    model = ADMETModel(**kw)
    smiles = list(dict.fromkeys(smiles))
    parts = []
    for i in range(0, len(smiles), chunk):
        log(f"  ADMET-AI {i + 1}-{min(i + chunk, len(smiles))}/{len(smiles)}")
        with open(os.devnull, "w") as dn, contextlib.redirect_stderr(dn), contextlib.redirect_stdout(io.StringIO()):
            parts.append(model.predict(smiles=smiles[i:i + chunk]))
    return pd.concat(parts) if parts else pd.DataFrame()
