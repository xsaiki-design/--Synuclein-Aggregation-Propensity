"""Extract Tables S3–S6 — the sequence signal — from the SI PDF.

S3 is the whole protein; S4/S5/S6 are the same 26 mutation types restricted to
the N-terminus, NAC and C-terminus. Together they say which *kinds* of
substitution change aggregation, and where in the sequence that holds.

Encoding note: these pages use the same CID-mapped numeric font as S7–S10 but
NOT always the same glyph offset — S5 and S6 are shifted by one, so a fixed map
decodes "1.000" as "20111". The offset is therefore calibrated per page by
trying candidates and keeping the one whose every cell decodes to a value
inside the known bounds (average in [-1,1], SD in [0,1], support a non-negative
integer). If no candidate validates, the page raises rather than guessing.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

RAW_PDF = Path("data/raw/SI.pdf")
OUT_CSV = Path("data/type_scores.csv")

_BASE = {1: " ", 2: "-", 3: ".", 4: "0", 5: "1", 6: "2", 7: "3",
         8: "4", 9: "5", 10: "6", 11: "7", 12: "8", 13: "9"}

# The 26 mutation-type rows, in the exact printed order shared by S3–S6.
TYPES = [
    "+ -> -", "- -> +", "nP -> P", "P -> nP", "Sm -> Lg", "Lg -> Sm",
    "nB -> B", "B -> nB", "nβ -> β", "β -> nβ",  # nβ/β rows
    "Pro -> X", "X -> Pro", "Gly -> X", "X -> Gly", "PG -> X", "X -> PG",
    "+ Cntrl", "- Cntrl", "nP Cntrl", "P Cntrl", "Sm Cntrl", "Lg Cntrl",
    "nB Cntrl", "B Cntrl", "nβ Cntrl", "β Cntrl",
]
SCOPES = {3: "whole_protein", 4: "N-terminus", 5: "NAC", 6: "C-terminus"}


def _decode(cell: str | None, shift: int) -> str:
    """Decode a cell. Glyphs outside the numeric map (letters, the trailing
    space) are dropped as separators; a wrong offset is caught by the bounds
    check and the spot checks below, not by these characters."""
    if not cell:
        return ""
    g = re.findall(r"\(cid:(\d+)\)", cell)
    if g:
        return "".join(c for c in (_BASE.get(int(x) - shift) for x in g) if c).strip()
    return cell.strip()


def _row_values(row, shift):
    """Read the four printed columns positionally: Type, Average, SD, Support.
    Scanning every cell for numbers would pick up digits out of the CID-mapped
    label column, so the column positions are used instead."""
    if len(row) < 4:
        return None
    avg, sd, sup = (_decode(row[1], shift), _decode(row[2], shift),
                    _decode(row[3], shift))
    return (avg, sd, sup)


_NUM = re.compile(r"-?\d+\.?\d*")


def _valid(parsed) -> bool:
    """Support is always an integer; average/SD are either both present and in
    range, or both blank (the SI leaves them empty when support is 0)."""
    for v in parsed:
        if v is None:
            return False
        avg, sd, sup = v
        if not sup.isdigit():
            return False
        if avg == "" and sd == "":
            if int(sup) != 0:
                return False
            continue
        if not (_NUM.fullmatch(avg) and _NUM.fullmatch(sd)):
            return False
        if not (-1.0 <= float(avg) <= 1.0 and 0.0 <= float(sd) <= 1.0):
            return False
    return True


def _parse_page(page) -> list[list[str]]:
    body = page.extract_tables()[0][1:]
    if len(body) != len(TYPES):
        raise ValueError(f"expected {len(TYPES)} rows, got {len(body)}")
    for shift in (0, 1, 2, -1):
        parsed = [_row_values(r, shift) for r in body]
        if _valid(parsed):
            return parsed
    raise ValueError("no CID offset decodes this page within known bounds")


# Values read off the rendered SI pages. The extraction must reproduce them
# exactly, or the CID offset has drifted and the result is not trustworthy.
_SPOT_CHECKS = [
    ("whole_protein", "+ -> -", -0.500, 0.866, 4),
    ("whole_protein", "X -> PG", -0.853, 0.354, 34),
    ("whole_protein", "B -> nB", -0.633, 0.706, 30),
    ("N-terminus", "X -> Pro", -0.813, 0.390, 16),
    ("N-terminus", "nB -> B", 0.667, 0.745, 6),
    ("NAC", "nP -> P", -0.800, 0.400, 10),
    ("NAC", "X -> PG", -0.889, 0.314, 18),
    ("C-terminus", "Pro -> X", 1.000, 0.000, 5),
    ("C-terminus", "B -> nB", -1.000, 0.000, 7),
]


def _verify(df: pd.DataFrame) -> None:
    for scope, typ, avg, sd, sup in _SPOT_CHECKS:
        r = df[(df.scope == scope) & (df["type"] == typ)]
        if len(r) != 1:
            raise ValueError(f"spot check missing: {scope}/{typ}")
        r = r.iloc[0]
        if not (abs(r.average - avg) < 1e-9 and abs(r.sd - sd) < 1e-9
                and int(r.support) == sup):
            raise ValueError(
                f"spot check failed for {scope}/{typ}: got "
                f"({r.average}, {r.sd}, {r.support}), expected ({avg}, {sd}, {sup})"
            )


def build(pdf_path: Path = RAW_PDF, out: Path = OUT_CSV) -> pd.DataFrame:
    import pdfplumber

    if not pdf_path.exists():
        raise FileNotFoundError(f"{pdf_path} not found; see README for the download step.")
    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        for idx, scope in SCOPES.items():
            for name, (avg, sd, sup) in zip(TYPES, _parse_page(pdf.pages[idx])):
                rows.append({
                    "scope": scope, "type": name,
                    "average": float(avg) if avg else None,
                    "sd": float(sd) if sd else None,
                    "support": int(sup),
                })
    df = pd.DataFrame(rows)
    _verify(df)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows) [REAL — Tables S3–S6, bounds-checked]")
    return df


if __name__ == "__main__":
    build()
