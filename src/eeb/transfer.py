"""Does a published correction transfer to data it was not fitted to?

A genome-scale model gets corrected: somebody finds reactions it is missing or gene associations it has
wrong, fixes them against a dataset, and publishes the repaired model. The repair is then reported as an
improvement, measured on the dataset it was made against. That is the only place it could have been
measured, and it is also the one place where the measurement means least.

This module asks the next question. The three essentiality experiments here were not used to make the
corrections in iML1515_corrected, which were derived from RB-TnSeq carbon-source fitness data. So grading
the corrected model against them is an out-of-sample test by construction, with no need to trust anybody's
discipline about held-out data.

Two things are needed to make the answer credible rather than merely negative:

  `changed_calls` — which genes actually changed, and which way. A difference of 0.03 in a summary
  statistic is not interpretable on its own; thirteen named genes in two pathways is.

  `bootstrap_delta` — whether a difference that small survives resampling the genes. Without this, a
  negative result is indistinguishable from noise, and reporting it would be as careless as reporting
  the positive one.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .gpr import describe_change, essential_genes
from .model import score

B_NUMBER = re.compile(r"b\d{4}")


def corrections_diff(before, after) -> dict:
    """What one model actually changes relative to another, as reactions, gene rules and bounds.

    A transfer result is only interpretable next to the content of the change, and a published model diff
    is rarely reported.

    Gene rules are compared by **Boolean equivalence**, not by which genes appear in them. The two are
    not the same test and the difference is the whole point: `A and B` and `A or B` contain identical
    gene sets and mean opposite things — a complex in which both genes are essential, against a pair of
    isozymes in which neither is. A gene-set comparison calls that pair unchanged. Since a rule change is
    only interesting here because of what it does to essentiality, every rule whose text differs is
    classified by its effect on the set of genes the rule makes essential:

      equivalent      the same Boolean function, written differently — bracketing, ordering
      more essential  the rule now makes strictly more genes essential (an isozyme was dropped)
      less essential  strictly fewer (an isozyme was added, or a complex relaxed)
      both            some genes gained essentiality and others lost it
      restructured    the function changed without changing which genes it makes essential
    """
    ra = {r.id: r for r in before.reactions}
    rb = {r.id: r for r in after.reactions}
    shared = set(ra) & set(rb)

    rules, equivalent_count = [], 0
    for rid in sorted(shared):
        old, new = ra[rid].gene_reaction_rule, rb[rid].gene_reaction_rule
        if old == new:
            continue
        kind = describe_change(old, new)
        if kind == "equivalent":
            equivalent_count += 1
            continue
        genes_old, genes_new = set(B_NUMBER.findall(old)), set(B_NUMBER.findall(new))
        rules.append({
            "reaction": rid, "change": kind,
            "dropped": " ".join(sorted(genes_old - genes_new)),
            "added": " ".join(sorted(genes_new - genes_old)),
            "newly_essential": " ".join(sorted(essential_genes(new) - essential_genes(old))),
            "no_longer_essential": " ".join(sorted(essential_genes(old) - essential_genes(new))),
            "rule_before": old, "rule_after": new,
        })
    return {
        "reactions_added": sorted(set(rb) - set(ra)),
        "reactions_removed": sorted(set(ra) - set(rb)),
        "gene_rules_changed": pd.DataFrame(rules),
        "gene_rules_equivalent": equivalent_count,
        "bounds_changed": sum(ra[r].bounds != rb[r].bounds for r in shared),
    }


def changed_calls(before: pd.Series, after: pd.Series, observed: pd.DataFrame,
                  sources: tuple[str, ...]) -> pd.DataFrame:
    """The genes whose essentiality call differs between two models, with what each experiment says.

    One row per changed gene, plus a `moved` column saying whether the change agrees with the majority of
    the experiments that have an opinion. This is the gene-level view of a summary statistic, and it is
    where the mechanism of a transfer failure becomes visible.
    """
    idx = before.index.intersection(after.index)
    differs = (before.reindex(idx) != after.reindex(idx)).fillna(False)
    changed = idx[differs]
    out = pd.DataFrame({"before": before.reindex(changed), "after": after.reindex(changed)})
    out = out.join(observed[["gene", *sources]])
    votes = observed.loc[changed, list(sources)]
    n_opinions = votes.notna().sum(axis=1)
    n_essential = votes.fillna(False).sum(axis=1)
    majority = (n_essential * 2 > n_opinions).where(n_opinions > 0)
    out["experiments_mostly_say"] = majority
    out["moved"] = np.where(out["after"] == majority, "toward",
                            np.where(out["before"] == majority, "away", "unclear"))
    return out


def direction(changed: pd.DataFrame, sources: tuple[str, ...]) -> pd.DataFrame:
    """For each experiment, how many of the changed calls became right and how many became wrong."""
    rows = []
    for src in sources:
        s = changed.dropna(subset=[src])
        rows.append({
            "experiment": src, "n_changed_with_opinion": len(s),
            "now_right": int(((s["after"] == s[src]) & (s["before"] != s[src])).sum()),
            "now_wrong": int(((s["before"] == s[src]) & (s["after"] != s[src])).sum()),
        })
    return pd.DataFrame(rows)


def bootstrap_delta(before: pd.Series, after: pd.Series, observed: pd.Series,
                    draws: int = 600, seed: int = 0) -> dict:
    """Resample genes to see whether a difference in MCC survives.

    Genes are the unit of resampling, which is the right unit here: the question is whether the corrected
    model would still have scored worse on another sample of genes, not on another sample of experiments.
    It is a crude interval — the genes are not independent, since they sit in pathways that the corrections
    act on together, so the real interval is wider than this one. It is reported because a bootstrap that
    *fails* to separate the two models settles the question anyway, and this one mostly does not.
    """
    rng = np.random.default_rng(seed)
    common = before.index.intersection(after.index).intersection(observed.dropna().index)
    deltas = np.empty(draws)
    for i in range(draws):
        draw = rng.choice(np.asarray(common), size=len(common), replace=True)
        deltas[i] = (score(after.reindex(draw), observed.reindex(draw))["mcc"]
                     - score(before.reindex(draw), observed.reindex(draw))["mcc"])
    return {
        "delta_mcc_mean": float(deltas.mean()),
        "delta_mcc_lo": float(np.quantile(deltas, 0.05)),
        "delta_mcc_hi": float(np.quantile(deltas, 0.95)),
        "p_worse": float((deltas < 0).mean()),
        "n_genes": len(common), "draws": draws,
    }


def reconcile(grades: pd.DataFrame, model_before: str, model_after: str) -> pd.DataFrame:
    """Why the raw error count and the MCC can disagree about the same change.

    Counting errors treats a lost true positive and a lost false positive as the same size of mistake.
    They are not. Essential genes are about a twelfth of the graded set, so a true positive is scarce and
    a false positive is cheap; a change that sheds more false alarms than it loses real calls reduces the
    error count while making the model worse at the only thing anybody wants it for. MCC prices the two
    classes by how common they are, which is why it moves the other way.

    Returns the trade explicitly, so the reader can check the arithmetic rather than take the claim.
    """
    rows = []
    for (medium, src), g in grades.groupby(["medium", "experiment"], sort=False):
        b = g[g["model"] == model_before]
        a = g[g["model"] == model_after]
        if b.empty or a.empty:
            continue
        b, a = b.iloc[0], a.iloc[0]
        rows.append({
            "medium": medium, "experiment": src,
            "prevalence": b["observed_essential"] / b["n"],
            "tp_lost": int(b["tp"] - a["tp"]), "fp_shed": int(b["fp"] - a["fp"]),
            "errors_before": int(b["fp"] + b["fn"]), "errors_after": int(a["fp"] + a["fn"]),
            "delta_errors": int((a["fp"] + a["fn"]) - (b["fp"] + b["fn"])),
            "delta_recall": float(a["recall"] - b["recall"]),
            "delta_precision": float(a["precision"] - b["precision"]),
            "delta_mcc": float(a["mcc"] - b["mcc"]),
        })
    return pd.DataFrame(rows)
