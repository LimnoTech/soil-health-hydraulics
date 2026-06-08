"""Fetch UNSODA 2.0 and extract a tidy bulk-density / organic-matter table.

UNSODA 2.0 (Nemes et al., 2001) is distributed as a Microsoft Access database
(`unsoda.mdb`) inside a ZIP on Ag Data Commons / figshare. This script downloads it,
exports the `soil_properties` and `general` tables with `mdbtools`, joins them, and
writes `data_temp/unsoda_bd_om.csv` — the subset of samples that report BOTH bulk
density and organic-matter content, with organic carbon computed via the van Bemmelen
factor (OC = 0.58 * OM).

`data_temp/` is git-ignored, so run this once to (re)generate the working data used by
the UNSODA cell in `2_organic_matter_water_holding.ipynb`.

Requires `mdbtools` on PATH (`mdb-export`). On macOS: `brew install mdbtools`.

    pixi run python notebooks/fetch_unsoda.py
"""
from __future__ import annotations

import io
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

UNSODA_ZIP_URL = "https://ndownloader.figshare.com/files/44334065"  # unsoda.ZIP on Ag Data Commons
VAN_BEMMELEN = 0.58  # OC = 0.58 * OM
DATA_DIR = Path(__file__).resolve().parent / "data_temp"


def main() -> None:
    if shutil.which("mdb-export") is None:
        sys.exit("mdb-export not found — install mdbtools (e.g. `brew install mdbtools`) and retry.")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    mdb_path = DATA_DIR / "unsoda.mdb"

    if not mdb_path.exists():
        print(f"Downloading UNSODA 2.0 from {UNSODA_ZIP_URL} ...")
        req = urllib.request.Request(UNSODA_ZIP_URL, headers={"User-Agent": "Mozilla/5.0"})
        blob = urllib.request.urlopen(req, timeout=180).read()
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            z.extractall(DATA_DIR)
        print(f"  extracted to {DATA_DIR}")

    def export(table: str) -> pd.DataFrame:
        out = subprocess.run(
            ["mdb-export", str(mdb_path), table], capture_output=True, text=True, check=True
        ).stdout
        return pd.read_csv(io.StringIO(out))

    sp = export("soil_properties")  # code, bulk_density, OM_content, ...
    gen = export("general")  # code, texture, horizon, depth_upper, depth_lower, ...

    both = sp[sp["bulk_density"].notna() & sp["OM_content"].notna()].merge(
        gen[["code", "texture", "horizon", "depth_upper", "depth_lower"]], on="code", how="left"
    )
    both = both.rename(columns={"bulk_density": "bulk_density_g_cm3"})
    both["OC_pct"] = VAN_BEMMELEN * both["OM_content"]
    both["is_mineral"] = both["OM_content"] <= 20  # exclude peat/organic outliers (OM > 20%)

    cols = [
        "code", "texture", "horizon", "depth_upper", "depth_lower",
        "bulk_density_g_cm3", "OM_content", "OC_pct", "is_mineral",
    ]
    out_path = DATA_DIR / "unsoda_bd_om.csv"
    both[cols].to_csv(out_path, index=False)
    print(f"Wrote {out_path}  ({len(both)} samples with both BD and OM)")


if __name__ == "__main__":
    main()
