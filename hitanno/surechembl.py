"""First patent per compound from the SureChEMBL bulk parquet files (remote or local) via DuckDB."""
import re

import pandas as pd

REMOTE = "https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/latest/"


def _clean_assignee(names):
    if names is None or len(names) == 0:
        return None
    return re.sub(r"\s*\([A-Z]{2}\)\s*$", "", str(names[0])).strip()


DEFAULT_CLASSES = r'^(C07|A61K|A61P)'


def first_patents(inchikeys, source=REMOTE, batch_size=5000, class_regex=DEFAULT_CLASSES, log=print):
    """Earliest-published patent per InChIKey.

    Exact InChIKey matches are preferred; if a compound only appears without matching
    stereo layers, the earliest patent for the same connectivity block is used and
    `patent_match` is 'connectivity'. SureChEMBL carries the publication date only.
    Text-mined name-to-structure errors produce unrelated old patents (e.g. firearms), so only
    patents with an IPC class matching `class_regex` are considered ('' disables the filter).
    """
    import duckdb
    src = source if source.endswith("/") else source + "/"
    con = duckdb.connect()
    if src.startswith("http"):
        con.execute("INSTALL httpfs; LOAD httpfs;")
    keys = sorted(set(inchikeys))
    frames = []
    allc = "list_concat(list_concat(coalesce(pt.ipc, []::VARCHAR[]), coalesce(pt.ipcr, []::VARCHAR[])), coalesce(pt.cpc, []::VARCHAR[]))"
    rx = class_regex.replace("'", "''")
    # keep a patent if any of its classes matches, or if it carries no classification at all
    cls = f"AND (len({allc}) = 0 OR len(list_filter({allc}, x -> regexp_matches(x, '{rx}'))) > 0)" if class_regex else ""
    for i in range(0, len(keys), batch_size):
        chunk = keys[i:i + batch_size]
        log(f"  SureChEMBL batch {i // batch_size + 1}/{-(-len(keys) // batch_size)} ({len(chunk)} compounds)")
        con.execute("CREATE OR REPLACE TEMP TABLE k AS SELECT * FROM (VALUES " +
                    ",".join(f"('{k}', '{k[:14]}')" for k in chunk) + ") t(ik, blk)")
        frames.append(con.execute(f"""
            WITH c AS (
                SELECT c.id AS cid, k.ik AS ik, (c.inchi_key = k.ik) AS exact
                FROM read_parquet('{src}compounds.parquet') c
                JOIN k ON left(c.inchi_key, 14) = k.blk),
            m AS (
                SELECT DISTINCT c.ik, c.exact, p.patent_id
                FROM read_parquet('{src}patent_compound_map.parquet') p JOIN c ON p.compound_id = c.cid)
            SELECT m.ik, m.exact, pt.patent_number, pt.publication_date, pt.assignee, pt.title, {allc} AS ipc
            FROM m JOIN read_parquet('{src}patents.parquet') pt ON pt.id = m.patent_id
            WHERE pt.publication_date IS NOT NULL
              {cls}
            QUALIFY row_number() OVER (PARTITION BY m.ik ORDER BY m.exact DESC, pt.publication_date, pt.patent_number) = 1
        """).df())
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    return pd.DataFrame({
        "inchikey": df["ik"],
        "patent_number": df["patent_number"],
        "patent_pub_date": pd.to_datetime(df["publication_date"]).dt.strftime("%Y-%m-%d"),
        "patent_originator": df["assignee"].map(_clean_assignee),
        "patent_use": df["title"],
        "patent_ipc": df["ipc"].map(lambda x: "; ".join(list(x)[:3]) if x is not None and len(x) else None),
        "patent_match": df["exact"].map({True: "exact", False: "connectivity"}),
    })
