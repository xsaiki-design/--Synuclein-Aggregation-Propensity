"""Extract the REAL condition tables (S7–S10) from the Pancoe et al. 2022 SI PDF.

What this SI actually contains (confirmed by reading all 12 pages):
  * The outcome is a CATEGORICAL aggregation score in {-1, 0, +1}
    (reduces / no effect / accelerates, normalised to the WT control) — NOT a
    continuous T½. T½ is not tabulated in this file.
  * Tables S3–S6 give scores averaged by mutation *type* and domain.
  * Tables S7–S10 give scores for the 6 familial mutants × condition category
    (buffer / salt / shaking speed / αS concentration), each with a support n.
    These are the only per-(variant × condition) data in the file, so they are
    the conditions-vs-sequence dataset.

Extraction detail: the label column is a CID-mapped font with no ToUnicode, so
the row *labels* are taken from the rendered tables (hard-coded below, in the
exact row order of each table). The *numbers* are decoded straight from the
PDF's own text layer via the font's fixed CID offset, and the row counts are
asserted to match. No value is invented; a mismatch raises.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

RAW_PDF = Path("data/raw/SI.pdf")
OUT_CSV = Path("data/conditions_familial.csv")

# Fixed CID -> character map for the numeric font used in Tables S7–S10
# (digit d is glyph d+4; '-'->2, '.'->3, ' '->1). Verified against every cell.
_CID = {1: " ", 2: "-", 3: ".", 4: "0", 5: "1", 6: "2", 7: "3",
        8: "4", 9: "5", 10: "6", 11: "7", 12: "8", 13: "9"}

FAMILIAL = ["A30P", "E46K", "H50Q", "G51D", "A53E", "A53T"]

# (table page index, dimension name, ordered (mutant, category) row labels).
# Row order is exactly as printed in the SI; counts are asserted at parse time.
_SPEC = {
    8:  ("buffer", [(m, c) for m in FAMILIAL for c in ("PBS", "Tris", "Other")]),
    9:  ("salt",   [(m, c) for m in FAMILIAL for c in ("High", "Low", "None")]),
    10: ("shaking", [(m, c) for m in FAMILIAL
                     for c in ("Fast", "Medium", "Slow", "None")]),
    # S10 (αS concentration) has an irregular order/last-row count in the SI:
    11: ("concentration",
         [("A30P", "High"), ("A30P", "Normal"), ("A30P", "Low"),
          ("E46K", "High"), ("E46K", "Normal"), ("E46K", "Low"),
          ("H50Q", "High"), ("H50Q", "Normal"), ("H50Q", "Low"),
          ("G51D", "High"), ("G51D", "Normal"), ("G51D", "Low"),
          ("A53E", "Normal"), ("A53E", "Low"), ("A53E", "High"),
          ("A53T", "Normal"), ("A53T", "Low")]),
}


def _decode(cell: str | None) -> str:
    if not cell:
        return ""
    glyphs = re.findall(r"\(cid:(\d+)\)", cell)
    if glyphs:
        return "".join(_CID.get(int(g), "?") for g in glyphs).strip()
    return cell.strip()  # already-plain cells (e.g. S3 page)


def _row_numbers(cells) -> list[str]:
    out = []
    for c in cells:
        d = _decode(c)
        if re.fullmatch(r"-?\d+\.?\d*", d):
            out.append(d)
    return out


# Values read off the rendered SI pages (S7-S10), one per table. The extraction
# must reproduce them or the CID offset has drifted; see ingest_types for the
# same guard on S3-S6.
_SPOT_CHECKS = [
    ("A30P", "buffer", "PBS", -0.75, 4),
    ("A53T", "buffer", "Tris", 0.94, 16),
    ("E46K", "salt", "None", 0.56, 9),
    ("G51D", "shaking", "Fast", -1.00, 4),
    ("A30P", "concentration", "High", -0.86, 7),
]


def _verify(df: pd.DataFrame) -> None:
    for mut, dim, cat, score, sup in _SPOT_CHECKS:
        r = df[(df.mutant == mut) & (df.dimension == dim) & (df.category == cat)]
        if len(r) != 1:
            raise ValueError(f"spot check missing: {mut}/{dim}/{cat}")
        r = r.iloc[0]
        if abs(r.score - score) > 1e-9 or int(r.support) != sup:
            raise ValueError(
                f"spot check failed for {mut}/{dim}/{cat}: got "
                f"({r.score}, {r.support}), expected ({score}, {sup})")


def build(pdf_path: Path = RAW_PDF, out: Path = OUT_CSV) -> pd.DataFrame:
    import pdfplumber

    if not pdf_path.exists():
        raise FileNotFoundError(
            f"{pdf_path} not found. Download the open-access SI from "
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159/ "
            "('Supporting_Information.pdf') and save it there."
        )
    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        for idx, (dim, labels) in _SPEC.items():
            table = pdf.pages[idx].extract_tables()[0][1:]  # drop header row
            nums = [_row_numbers(r) for r in table]
            nums = [n for n in nums if len(n) == 2]  # (average, support)
            if len(nums) != len(labels):
                raise ValueError(
                    f"{dim}: parsed {len(nums)} numeric rows but expected "
                    f"{len(labels)} labels — table layout changed, refusing to guess."
                )
            for (mut, cat), (avg, sup) in zip(labels, nums):
                rows.append({"mutant": mut, "dimension": dim, "category": cat,
                             "score": float(avg), "support": int(sup)})
    df = pd.DataFrame(rows)
    _verify(df)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows) [REAL — Tables S7–S10, verified]")
    return df


if __name__ == "__main__":
    build()
