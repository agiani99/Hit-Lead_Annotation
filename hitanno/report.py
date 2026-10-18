"""Radar plots (green = allowed area) and the static HTML report."""
import base64
import html
import io
import datetime
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

MEAS, PRED, GREEN = "#111111", "#1565c0", "#2e7d32"
R_THR = 0.6  # every axis threshold is drawn at this radius, so the allowed area is a disc


def _scale(ax, source):
    """(direction, thr, lo, hi) on the scale of the value's source."""
    if source == "predicted" and ax["kind"] != "admet":
        return ax["direction"], ax["pred_thr"], ax["pred_lo"], ax["pred_hi"]
    return ax["direction"], ax.get("thr"), ax["lo"], ax["hi"]


def radius(ax, v, source):
    d, thr, lo, hi = _scale(ax, source)
    if d == "band":
        mid, half = (ax["lo_ok"] + ax["hi_ok"]) / 2, (ax["hi_ok"] - ax["lo_ok"]) / 2
        dist = abs(v - mid)
        r = R_THR * dist / half if dist <= half else R_THR + (1 - R_THR) * (dist - half) / half
    else:
        if d == "higher_better":
            v, thr, lo, hi = -v, -thr, -hi, -lo
        r = R_THR * (v - lo) / (thr - lo) if v <= thr else R_THR + (1 - R_THR) * (v - thr) / (hi - thr)
    return float(min(max(r, 0.0), 1.0))


def radar_png(row, axes, title, originator, date, expired):
    n = len(axes)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    fig = plt.figure(figsize=(4.4, 4.6), dpi=110)
    ax = fig.add_axes([0.14, 0.08, 0.72, 0.68], projection="polar")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["polar"].set_color("#bbbbbb")
    ax.grid(color="#dddddd", lw=0.6)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.fill(th, np.full_like(th, R_THR), color="#43a047", alpha=0.22, lw=0)
    ax.plot(th, np.full_like(th, R_THR), color="#2e7d32", lw=1)

    pts = []
    for a, t in zip(ang, axes):
        v, s = row.get(t["key"]), row.get(f"{t['key']}_source", "")
        pts.append((a, radius(t, v, s), s) if v is not None and not pd.isna(v) else (a, None, ""))
    have = [p for p in pts if p[1] is not None]
    if len(have) > 2:
        xs = [p[0] for p in have] + [have[0][0]]
        ys = [p[1] for p in have] + [have[0][1]]
        ax.plot(xs, ys, color="#555555", lw=1.2)
        ax.fill(xs, ys, color="#555555", alpha=0.08)
    for a, r, s in have:
        c = PRED if s == "predicted" else MEAS
        ax.plot([a], [r], "o", ms=5.5, color=c, zorder=5)

    ax.set_xticks(ang)
    ax.set_xticklabels([])
    for a, t, (_, r, _s) in zip(ang, axes, pts):
        ax.text(a, 1.16, t["label"], ha="center", va="center", fontsize=7.5,
                color="#222222" if r is not None else "#aaaaaa")
        if r is not None and r > R_THR:
            ax.text(a, 1.03, "!", ha="center", va="center", fontsize=8, color="#c62828", fontweight="bold")
    fig.text(0.5, 0.955, title, ha="center", va="center", fontsize=10.5, fontweight="bold", color="#111111")
    fig.text(0.5, 0.913, originator, ha="center", va="center", fontsize=8, color="#555555")
    fig.text(0.5, 0.875, date, ha="center", va="center", fontsize=8,
             color=GREEN if expired else "#555555", fontweight="bold" if expired else "normal")
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _fmt(v, nd=2):
    return "n/a" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{nd}f}"


def build_report(df, cfg, out_html, min_phase=1, cols=3, patent_term=25, log=print):
    axes = [a for a in cfg if a["key"] in df and df[a["key"]].notna().any()]
    dropped = [a["label"] for a in cfg if a not in axes]
    if dropped:
        log(f"  axes with no data for any compound (dropped): {', '.join(dropped)}")
    sel = df[df.max_phase >= min_phase].sort_values(["max_phase", "pchembl"], ascending=[False, False])
    sections = []
    for phase, grp in sel.groupby("max_phase", sort=False):
        cards = []
        for r in grp.to_dict("records"):
            org = r["patent_originator"] if isinstance(r["patent_originator"], str) else "originator n/a"
            date, expired = "date n/a", False
            if isinstance(r["patent_pub_date"], str):
                d = datetime.date.fromisoformat(r["patent_pub_date"])
                expired = d.replace(year=d.year + patent_term) < datetime.date.today()
                date = d.isoformat() + (f"  (expired, >{patent_term}y)" if expired else "")
            png = radar_png(r, axes, str(r["name"])[:38], org[:60], date, expired)
            n_pred = sum(r.get(f"{a['key']}_source") == "predicted" for a in axes)
            cap = (f"{html.escape(str(r['chembl_id']))} · pChEMBL {_fmt(r['pchembl'])} · {n_pred} predicted axes"
                   + (f"<br><span class='use'>{html.escape(str(r['patent_use']))[:110]}</span>"
                      if r.get("patent_use") and not pd.isna(r["patent_use"]) else ""))
            cards.append(f"<figure><img src='data:image/png;base64,{png}'><figcaption>{cap}</figcaption></figure>")
        sections.append(f"<h2>Phase {int(phase)} <small>{len(cards)} compound(s)</small></h2>"
                        f"<div class='grid'>{''.join(cards)}</div>")
    doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>HIT annotation radar plots</title><style>
body{{font-family:system-ui,Segoe UI,sans-serif;margin:24px auto;max-width:1400px;color:#222;padding:0 16px}}
h2{{border-bottom:2px solid #ddd;padding-bottom:4px;margin-top:36px}} h2 small{{font-weight:400;color:#777;font-size:.6em}}
.grid{{display:grid;grid-template-columns:repeat({cols},1fr);gap:16px}}
figure{{margin:0;border:1px solid #e3e3e3;border-radius:8px;padding:8px;text-align:center}}
figure img{{max-width:100%}} figcaption{{font-size:12px;color:#555}} .use{{color:#888}}
.legend{{background:#f6f8f6;border-radius:8px;padding:10px 14px;font-size:13px}}
@media(max-width:900px){{.grid{{grid-template-columns:1fr}}}}</style></head><body>
<h1>HIT annotation radar plots</h1>
<div class='legend'>Centre = good, rim = bad. <b style='color:#2e7d32'>Green disc = allowed area</b>: a point inside it passes
the threshold, a red <b style='color:#c62828'>!</b> marks a violation. <b>Black dot = measured / ChEMBL value</b>,
<b style='color:#1565c0'>blue dot = predicted</b>; grey axis label = no data. Strictest value is used when several exist.
Dates are first patent <i>publication</i> dates (SureChEMBL); a <b style='color:#2e7d32'>green date</b> means date + {patent_term} years is already past (patent term expired).</div>
{''.join(sections) or '<p>No compounds at or above the minimum phase.</p>'}</body></html>"""
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(doc)
