"""Five published opinions about which E. coli genes are essential, put on the same index.

A genome-scale metabolic model's essentiality accuracy is always quoted against one experiment. There is
more than one experiment, they were done with different technologies by different labs, and they disagree.
This module loads them, maps them all to b-numbers, and makes the disagreement visible before any model is
asked to predict anything.

The sources, all of E. coli K-12 and four of them of the same strain, BW25113:

  keio      Baba et al. (2006). Systematic single-gene deletion. A gene is called essential when no
            verified deletion mutant could be made. Supplementary Table 3, column "verified" == 0.
  tradis    Goodall et al. (2018). Dense transposon mutagenesis, essentiality inferred from the gap in
            insertions. mBio Table S4, which also reports an "unclear" class — the only source honest
            enough to have one.
  pec       Profiling of the E. coli Chromosome. A curated database call, essential / non-essential /
            unknown.
  genetic   The classical genetic literature's call, as collated in the Keio supplement.
  tn5       An earlier MG1655 Tn5 transposon screen, likewise collated in the Keio supplement.

Technologies fail in different ways, and the ways matter. A transposon screen miscalls a gene whose
essential part is short, because insertions elsewhere in it are tolerated, and miscalls the gene upstream
of an essential one, because the insertion knocks out the promoter. A deletion collection miscalls a gene
whose mutant grows only with a suppressor mutation, and inherits any error in library construction. So
where two sources disagree, neither is automatically the one that is wrong, and a model graded against
either is being graded partly on that source's failure modes.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

RAW = Path(__file__).resolve().parents[2] / "data" / "raw"
SOURCES = ("keio", "tradis", "pec", "genetic", "tn5")
# The four that are calls about the same strain by comparable methods. `genetic` is a literature
# collation rather than one experiment, so it is loaded but kept out of the headline comparison.
EXPERIMENTAL = ("keio", "tradis", "pec", "tn5")


def _keio_table() -> pd.DataFrame:
    """Baba et al. Supplementary Table 3: one row per gene, with several sources' calls side by side."""
    raw = pd.read_excel(RAW / "keio_baba2006_table3.xls", sheet_name="Sup Table 3", header=None)
    body = raw.iloc[4:].reset_index(drop=True)
    out = pd.DataFrame({
        "gene": body[1].astype(str).str.strip(),
        "b_number": body[8].astype(str).str.strip().str.lower(),
        "verified": pd.to_numeric(body[14], errors="coerce"),
        "genetic_call": body[20].astype(str).str.strip(),
        "pec_call": body[21].astype(str).str.strip(),
        "tn5_call": body[22].astype(str).str.strip(),
        "keio_score": pd.to_numeric(body[23], errors="coerce"),
    })
    return out[out["b_number"].str.match(r"^b\d{4}$", na=False)].copy()


def _call(series: pd.Series, essential: set[str], non_essential: set[str]) -> pd.Series:
    """Map a source's labels onto True / False / NA, where NA means the source declined to say."""
    return series.map(lambda v: True if v in essential else (False if v in non_essential else pd.NA))


def _pec_table() -> pd.DataFrame:
    """The PEC export. Tab separated inside an RTF wrapper, so the wrapper is stripped first."""
    text = (RAW / "pec_genes.txt.rtf").read_text(encoding="latin-1", errors="ignore")
    text = re.sub(r"\\'[0-9a-f]{2}", "", text)
    text = re.sub(r"\{\\\*?[^{}]*\}|[{}]|\\[A-Za-z]+\d*\s?", " ", text)
    rows = []
    for line in text.split("\n"):
        parts = line.split("\t")
        if len(parts) < 4 or not parts[2].strip().isdigit():
            continue
        klass = int(parts[2].strip())
        b = next((a.strip().lower() for a in parts[3].split(",") if re.match(r"^b\d{4}$", a.strip())), None)
        if b is None or klass not in (1, 2, 3):
            continue
        rows.append({"b_number": b, "gene": parts[1].strip(),
                     "pec_file": {1: True, 2: False, 3: pd.NA}[klass]})
    return pd.DataFrame(rows).drop_duplicates("b_number")


def _tradis_table(keio: pd.DataFrame) -> pd.DataFrame:
    """Goodall et al. Table S4. Keyed by gene name, so it is mapped to b-numbers through the Keio table."""
    raw = pd.read_excel(RAW / "tradis_goodall2018_tableS4.xlsx", header=1)
    raw.columns = [str(c).strip() for c in raw.columns]
    raw["gene"] = raw["Gene"].astype(str).str.strip()
    essential = raw["Essential"].astype(str).str.lower() == "true"
    unclear = raw["Unclear"].astype(str).str.lower() == "true"
    raw["tradis"] = pd.Series([pd.NA if u else bool(e) for e, u in zip(essential, unclear)],
                              index=raw.index, dtype="boolean")
    names = keio.drop_duplicates("gene").set_index("gene")["b_number"]
    raw["b_number"] = raw["gene"].map(names)
    return (raw.dropna(subset=["b_number"])
            [["b_number", "tradis", "Insertion Index Score", "Log Likelihood Ratio"]]
            .rename(columns={"Insertion Index Score": "tradis_insertion_index",
                             "Log Likelihood Ratio": "tradis_llr"})
            .drop_duplicates("b_number"))


def load() -> pd.DataFrame:
    """One row per gene, one boolean column per source, indexed by b-number.

    True means that source called the gene essential, False that it called it non-essential, and NA that
    it had no opinion — which is a third state and is kept as one rather than folded into False.
    """
    keio = _keio_table()
    table = pd.DataFrame({
        "gene": keio["gene"].values,
        "b_number": keio["b_number"].values,
        "keio": (keio["verified"] == 0).where(keio["verified"].notna(), pd.NA),
        "genetic": _call(keio["genetic_call"], {"E"}, {"N", "N*"}).values,
        "pec": _call(keio["pec_call"], {"E"}, {"N"}).values,
        "tn5": _call(keio["tn5_call"], {"E"}, {"N"}).values,
        "keio_score": keio["keio_score"].values,
    }).drop_duplicates("b_number").set_index("b_number")

    table = table.join(_tradis_table(keio).set_index("b_number"))
    pec_file = _pec_table().set_index("b_number")["pec_file"]
    # the curated database export agrees with the collated column where both speak; it fills the gaps
    table["pec"] = table["pec"].fillna(pec_file.reindex(table.index))
    for source in SOURCES:
        table[source] = table[source].astype("boolean")
    return table


def agreement(table: pd.DataFrame, sources: tuple[str, ...] = EXPERIMENTAL) -> pd.DataFrame:
    """How often each pair of sources agrees, among the genes where both have an opinion."""
    rows = []
    for i, a in enumerate(sources):
        for b in sources[i + 1:]:
            both = table[[a, b]].dropna()
            agree = (both[a] == both[b]).sum()
            rows.append({
                "source_a": a, "source_b": b, "n_compared": len(both),
                "agree": int(agree), "agreement": agree / len(both) if len(both) else float("nan"),
                "essential_in_a_only": int((both[a] & ~both[b]).sum()),
                "essential_in_b_only": int((~both[a] & both[b]).sum()),
            })
    return pd.DataFrame(rows)


def consensus(table: pd.DataFrame, sources: tuple[str, ...] = EXPERIMENTAL) -> pd.DataFrame:
    """Per gene: how many sources call it essential, how many have an opinion, and whether they are
    unanimous. The unanimous genes are the only ones on which a model can be graded without the grade
    depending on which experiment you believe."""
    opinions = table[list(sources)]
    n_votes = opinions.notna().sum(axis=1)
    n_essential = opinions.fillna(False).sum(axis=1)
    out = pd.DataFrame({
        "gene": table["gene"], "n_sources": n_votes, "n_essential": n_essential,
        "n_non_essential": n_votes - n_essential,
    })
    out["unanimous"] = (n_votes > 1) & ((n_essential == 0) | (n_essential == n_votes))
    out["contested"] = (n_votes > 1) & ~out["unanimous"]
    out["consensus_essential"] = out["unanimous"] & (n_essential > 0)
    return out
