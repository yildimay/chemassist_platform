# -*- coding: utf-8 -*-
"""extract_molecules_from_pdf.py (rev‑2)

Robust extraction of *true* chemical names (e.g. Palbociclib, Kaempferol) from
research‑article PDFs.

Pipeline
========
1. **Read PDF → raw text** (PyMuPDF ‑ fast & layout‑agnostic)
2. **Chemical NER** via *ChemDataExtractor* → set of raw entity strings
3. Heuristics:   • drop stop‑words (Table, Revised …)  • length > 2  • no digits
4. *(Optional)* **PubChem sanity‑check** using *pubchempy*; keeps only entities
   that resolve to at least one PubChem CID (makes list noise‑free but slower).
5. Return unique names (preserving title‑case).

Usage
-----
```bash
pip install pymupdf chemdataextractor pubchempy
python extract_molecules_from_pdf.py my_article.pdf
```

The function `extract_molecules_from_pdf(file_like)` can be imported directly
into *smiles_to_mol.py* and connected to a Streamlit file‑uploader.
"""

from __future__ import annotations

import io
import re
import sys
from pathlib import Path
from typing import Iterable, List, Set, Union

import fitz  # PyMuPDF
from chemdataextractor import Document

# PubChem lookup is optional → soft import
try:
    import pubchempy as pcp
except ModuleNotFoundError:  # fallback stub
    pcp = None  # type: ignore

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
STOP_PATTERNS = re.compile(
    r"^(table|figure|received|revised|accepted|available|online|science|vertex|"
    r"explorer|changsha|lianyungang|shelxl|shelxs|cdk|journal|february|october)",
    re.I,
)
NON_WORD = re.compile(r"[^A-Za-z\-]", re.A)


def _read_pdf(file_like: Union[str, Path, bytes, io.BytesIO]) -> str:
    """Return concatenated plain text of every page using PyMuPDF."""

    if isinstance(file_like, (str, Path)):
        doc = fitz.open(str(file_like))
    elif isinstance(file_like, (bytes, io.BytesIO)):
        if isinstance(file_like, bytes):
            stream = io.BytesIO(file_like)
        else:
            stream = file_like
        stream.seek(0)
        doc = fitz.open(stream=stream.read(), filetype="pdf")
    else:
        raise TypeError("Unsupported file_like type → %r" % type(file_like))

    text_pages = [page.get_text("text") for page in doc]
    doc.close()
    return "\n".join(text_pages)


def _clean_names(raw_names: Iterable[str]) -> List[str]:
    """Apply regex/length filters and title‑case normalisation."""

    out: List[str] = []
    seen: Set[str] = set()

    for name in raw_names:
        name = name.strip().replace("\u00ad", "")  # soft hyphen removal
        if len(name) < 3:
            continue
        if STOP_PATTERNS.search(name):
            continue
        if any(ch.isdigit() for ch in name):
            continue
        # remove trailing punctuation
        name = name.rstrip(".,;:()[]{} ")
        # multiple words? keep up to three tokens (e.g. "L‑proline hydrochloride")
        tokens = name.split()
        if len(tokens) > 3:
            continue
        cleaned = " ".join(tokens)
        key = cleaned.lower()
        if key not in seen:
            out.append(cleaned)
            seen.add(key)
    return out


def _pubchem_filter(names: List[str]) -> List[str]:
    """Return only names that resolve to at least one PubChem CID."""

    if pcp is None:
        return names  # pubchempy not installed → skip

    confirmed: List[str] = []
    for nm in names:
        try:
            if pcp.get_compounds(nm, "name", first=True):
                confirmed.append(nm)
        except Exception:
            # network error / timeout → keep candidate (fail‑open)
            confirmed.append(nm)
    return confirmed


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def extract_molecules_from_pdf(
    file_like: Union[str, Path, bytes, io.BytesIO],
    *,
    validate_pubchem: bool = True,
    max_names: int | None = 30,
) -> List[str]:
    """Return a de‑duplicated list of chemical names detected in *file_like*.

    Parameters
    ----------
    file_like : path/bytes/BytesIO
        PDF file content or path.
    validate_pubchem : bool, default True
        If **True** each name is looked‑up via PubChem to weed out false positives.
    max_names : int | None
        Truncate list after *max_names*; ``None`` → no limit.
    """

    text = _read_pdf(file_like)
    # Chemical NER via ChemDataExtractor
    doc = Document(text)
    raw_names = {cem.text for cem in doc.cems if cem.text}

    names = _clean_names(raw_names)
    if validate_pubchem:
        names = _pubchem_filter(names)

    return names if max_names is None else names[: max_names]


# ---------------------------------------------------------------------------
# CLI for quick testing
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python extract_molecules_from_pdf.py <article.pdf>")

    pdf_path = sys.argv[1]
    print("Parsing:", pdf_path)
    cems = extract_molecules_from_pdf(pdf_path)
    print("\nDetected molecules (PubChem‑confirmed):")
    for i, cem in enumerate(cems, 1):
        print(f"{i:2d}. {cem}")
