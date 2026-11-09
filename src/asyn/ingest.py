"""Build data/aggregation_curated.csv from the real supplementary PDF.

Source (open access, author-sanctioned for computational reuse):
  Pancoe SX, Wang YJ, Shimogawa M, Perez RM, Giannakoulias S, Petersson EJ.
  "Effects of Mutations and Post-Translational Modifications on alpha-Synuclein
  In Vitro Aggregation." J Mol Biol 2022; 434(23):167859.
  PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159
  The tabulated outcome is T_half: time to half-maximal ThT fluorescence.

This module does NOT invent data. It extracts the tables from the PDF, dumps
them for human review, and maps the reviewed columns into the project schema.
Because the exact table layout can only be confirmed against the actual file,
the column mapping in `_map_rows` is a REVIEW POINT: run `--inspect` first,
confirm the dumped headers, then wire them up. Until SI.pdf is present and the
mapping is confirmed, this raises rather than guessing.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .data import REQUIRED_COLUMNS

RAW_PDF = Path("data/raw/SI.pdf")
DUMP_DIR = Path("data/raw")
CURATED_CSV = Path("data/aggregation_curated.csv")

DOWNLOAD_HINT = (
    f"{RAW_PDF} not found.\n"
    "Download the open-access supplementary PDF (one click, no login) from:\n"
    "  https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159/\n"
    "  -> 'Supporting_Information.pdf' (NIHMS1869546-supplement)\n"
    f"and save it as {RAW_PDF}. Then re-run `make data`."
)


def extract_tables(pdf_path: Path = RAW_PDF) -> list[pd.DataFrame]:
    """Return every table pdfplumber finds, one DataFrame per table."""
    import pdfplumber  # imported lazily so synth-only runs need no PDF libs

    if not pdf_path.exists():
        raise FileNotFoundError(DOWNLOAD_HINT)

    tables: list[pd.DataFrame] = []
    with pdfplumber.open(pdf_path) as pdf:
        for pno, page in enumerate(pdf.pages, start=1):
            for tno, tbl in enumerate(page.extract_tables(), start=1):
                if not tbl or len(tbl) < 2:
                    continue
                header, *body = tbl
                cols = [(_clean(c) or f"col{i}") for i, c in enumerate(header)]
                df = pd.DataFrame(body, columns=cols)
                df.insert(0, "_table", f"p{pno}_t{tno}")
                tables.append(df)
    return tables


def inspect(pdf_path: Path = RAW_PDF) -> Path:
    """Dump all extracted tables to one CSV for human review of the layout."""
    tables = extract_tables(pdf_path)
    if not tables:
        raise RuntimeError(
            f"No tables extracted from {pdf_path}. The PDF may be image-based; "
            "it then needs OCR, or copy the tables out by hand."
        )
    DUMP_DIR.mkdir(parents=True, exist_ok=True)
    out = DUMP_DIR / "si_tables_dump.csv"
    pd.concat(tables, ignore_index=True).to_csv(out, index=False)
    print(f"Extracted {len(tables)} tables -> {out}")
    print("Review the headers there, then confirm the mapping in _map_rows().")
    return out


def _map_rows(tables: list[pd.DataFrame]) -> pd.DataFrame:
    """Map reviewed PDF tables into the project schema.

    REVIEW POINT: fill this in against the confirmed `si_tables_dump.csv`
    headers. It must never fabricate a value; a condition the paper does not
    report for a row stays missing (NaN), which is itself recorded data.
    """
    raise NotImplementedError(
        "Column mapping not yet confirmed against the real PDF.\n"
        "Run `python -m asyn.ingest --inspect`, confirm the dumped table "
        "headers, then implement _map_rows() to emit REQUIRED_COLUMNS:\n  "
        + ", ".join(REQUIRED_COLUMNS)
    )


def build(pdf_path: Path = RAW_PDF, out: Path = CURATED_CSV) -> Path:
    tables = extract_tables(pdf_path)
    df = _map_rows(tables)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"mapped frame missing columns: {missing}")
    df = df[REQUIRED_COLUMNS]
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows) [REAL DATA from {pdf_path.name}]")
    return out


def _clean(s) -> str:
    return " ".join(str(s or "").split())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inspect", action="store_true",
                    help="extract and dump tables for review, don't build")
    args = ap.parse_args()
    if args.inspect:
        inspect()
    else:
        build()
