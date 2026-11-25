"""Grading a genome-scale metabolic model against an essentiality experiment.

The model predicts that deleting a gene is lethal when, with the reactions that gene encodes switched off,
the linear programme can no longer carry flux through the biomass reaction. That prediction depends on two
things besides the model: a growth threshold, and a medium.

**The medium is not a detail.** A gene for histidine biosynthesis is essential in minimal medium and
dispensable in rich medium, because in rich medium the cell can import histidine instead. iML1515's default
is M9 with glucose. The Keio collection was built in LB. Grading a minimal-medium model against rich-medium
essentiality therefore scores every amino-acid biosynthesis gene as a false positive, and the resulting
accuracy says more about the mismatch than about the model. Both media are therefore run, and the gap
between them is reported rather than a choice being made quietly.

**The models do not arrive on the same medium.** iML1515_corrected, as published, has nearly all of its
uptake reactions closed and is infeasible out of the box: the authors saved it part-way through an
analysis that set the medium itself. Grading the five models therefore starts by putting them all on the
same M9 glucose base, plus whatever uptakes each model's author deliberately left open — for the corrected
model that is biotin, pantothenate and thiamine, which is part of the correction rather than an accident.

**The threshold is not a detail either.** Deleting a gene rarely takes growth to exactly zero; it takes it
to some fraction of wild type. Where the cut is placed moves the essential set, so it is swept.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

MODELS = Path(__file__).resolve().parents[2] / "data" / "external" / "models"
MODEL_FILES = {
    "iJR904": "iJR904.xml", "iAF1260": "iAF1260.xml", "iJO1366": "iJO1366.xml",
    "iML1515": "iML1515.xml",
    "iML1515_corrected": "iML1515_model_adjusted_all_corrections.xml",
}

# The components of a rich medium, as exchange reaction ids. Twenty amino acids, four nucleosides and the
# vitamins and cofactors a rich broth supplies. Opening these is what makes a biosynthesis gene
# dispensable, which is the condition the Keio collection was built in.
RICH_SUPPLEMENTS = [
    "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e", "EX_asp__L_e", "EX_cys__L_e", "EX_gln__L_e",
    "EX_glu__L_e", "EX_gly_e", "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_pro__L_e", "EX_ser__L_e", "EX_thr__L_e", "EX_trp__L_e",
    "EX_tyr__L_e", "EX_val__L_e",
    "EX_adn_e", "EX_gsn_e", "EX_cytd_e", "EX_uri_e", "EX_thymd_e", "EX_ins_e",
    "EX_thm_e", "EX_pnto__R_e", "EX_nac_e", "EX_pydx_e", "EX_ribflv_e", "EX_btn_e", "EX_fol_e",
    "EX_4abz_e", "EX_cbl1_e", "EX_lipoate_e", "EX_hxan_e", "EX_ura_e",
]


def load_model(name: str):
    """Read one of the five SBML models. Slow (a few seconds), so callers should cache."""
    import cobra
    return cobra.io.read_sbml_model(str(MODELS / MODEL_FILES[name]))


# M9 with glucose, taken from iML1515's own default so that every model is compared on one base.
M9_GLUCOSE = {
    "EX_glc__D_e": -10.0, "EX_o2_e": -1000.0, "EX_nh4_e": -1000.0, "EX_pi_e": -1000.0,
    "EX_so4_e": -1000.0, "EX_h2o_e": -1000.0, "EX_h_e": -1000.0, "EX_co2_e": -1000.0,
    "EX_k_e": -1000.0, "EX_na1_e": -1000.0, "EX_cl_e": -1000.0, "EX_ca2_e": -1000.0,
    "EX_mg2_e": -1000.0, "EX_mn2_e": -1000.0, "EX_fe2_e": -1000.0, "EX_fe3_e": -1000.0,
    "EX_zn2_e": -1000.0, "EX_cu2_e": -1000.0, "EX_cobalt2_e": -1000.0, "EX_ni2_e": -1000.0,
    "EX_mobd_e": -1000.0, "EX_sel_e": -1000.0, "EX_slnt_e": -1000.0, "EX_tungs_e": -1000.0,
}


def set_medium(model, rich: bool, uptake: float = 10.0) -> dict:
    """Put the model on M9 glucose, keeping any extra uptake its author left open, and optionally
    supplement it to something like a rich broth.

    Returns the uptakes actually applied, since not every model generation has an exchange reaction for
    every component.
    """
    authors_extras = {r.id: r.lower_bound for r in model.exchanges
                      if r.lower_bound < 0 and r.id not in M9_GLUCOSE}
    for exchange in model.exchanges:
        exchange.lower_bound = 0.0
    applied = {}
    for ex, bound in {**M9_GLUCOSE, **authors_extras}.items():
        if ex in model.reactions:
            model.reactions.get_by_id(ex).lower_bound = bound
            applied[ex] = bound
    if rich:
        for ex in RICH_SUPPLEMENTS:
            if ex in model.reactions:
                model.reactions.get_by_id(ex).lower_bound = -uptake
                applied[ex] = -uptake
    return applied


def predict_essential(model, threshold: float = 0.01, processes: int = 1) -> pd.DataFrame:
    """Single-gene deletions. Returns one row per gene with the growth rate left after deleting it, as a
    fraction of wild type, and the resulting essentiality call."""
    from cobra.flux_analysis import single_gene_deletion
    wild_type = model.slim_optimize()
    results = single_gene_deletion(model, model.genes, processes=processes)
    genes = [next(iter(ids)) for ids in results["ids"]]
    growth = results["growth"].to_numpy()
    # An infeasible deletion comes back as NaN, and a lethal one often comes back as a tiny negative
    # number instead of exactly zero — LP solver tolerance, not biology. Both mean no growth.
    growth = np.clip(np.nan_to_num(growth, nan=0.0), 0.0, None)
    out = pd.DataFrame({
        "b_number": [str(g).strip().lower() for g in genes],
        "growth": growth,
        "growth_fraction": growth / wild_type if wild_type > 0 else np.nan,
    })
    if not (wild_type and wild_type > 1e-9):
        raise ValueError("the model does not grow on this medium, so nothing can be called essential")
    out["predicted_essential"] = out["growth_fraction"] < threshold
    return out.drop_duplicates("b_number").set_index("b_number")


def score(predicted: pd.Series, observed: pd.Series) -> dict:
    """Confusion matrix and the summary measures, over the genes where both the model and the experiment
    have something to say.

    Accuracy is reported because the literature reports it, and is nearly meaningless here: essential genes
    are a twentieth of the genome, so calling everything dispensable scores about 95%. Matthews correlation
    and the F1 over the essential class are what actually move.
    """
    both = pd.DataFrame({"p": predicted, "o": observed}).dropna()
    p = both["p"].astype(bool).to_numpy()
    o = both["o"].astype(bool).to_numpy()
    tp, fp = int((p & o).sum()), int((p & ~o).sum())
    tn, fn = int((~p & ~o).sum()), int((~p & o).sum())
    denom = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    return {
        "n": len(both), "observed_essential": int(o.sum()), "predicted_essential": int(p.sum()),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "accuracy": (tp + tn) / len(both) if len(both) else float("nan"),
        "precision": precision, "recall": recall,
        "f1": 2 * precision * recall / (precision + recall) if precision and recall else float("nan"),
        "mcc": (tp * tn - fp * fn) / denom if denom else float("nan"),
    }
