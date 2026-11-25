"""Grade the models against every essentiality experiment, and decompose where the disagreement goes.

  results/agreement.csv     how often each pair of experiments agrees
  results/consensus.csv     per gene: how many experiments call it essential, and whether they are unanimous
  results/grades.csv        every model x medium x experiment, with the confusion matrix
  results/split.csv         the same, split into genes the experiments agree on and genes they contest
  results/threshold.csv     how the grade moves with the growth threshold
  results/corrections.csv   the published corrections, graded on experiments they were not fitted to
  results/reconcile.csv     why the error count and the MCC disagree about those corrections
  results/changed.csv       the genes whose call the corrections changed, and which way
  results/bootstrap.csv     whether the change in MCC survives resampling the genes
  results/corrections_diff.csv  what the corrections actually change in the model
  results/common_genes.csv  every model graded on the genes all four of them contain
  results/strict_uptakes.csv  the same, with every model on strictly identical uptakes

Usage:  PYTHONPATH=src python scripts/02_grade.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eeb.essentiality import agreement, consensus, load
from eeb.model import MODEL_FILES, load_model, predict_essential, score, set_medium
from eeb.transfer import bootstrap_delta, changed_calls, corrections_diff, direction, reconcile

RESULTS = ROOT / "results"
SOURCES = ("keio", "tradis", "pec")
MEDIA = {"minimal": False, "rich": True}
BEFORE, AFTER = "iML1515", "iML1515_corrected"


def predictions(names=tuple(MODEL_FILES), strict: bool = False) -> dict:
    """Deletion scans for every model under both media. The slow part, and still only a minute or two."""
    out = {}
    for name in names:
        model = load_model(name)
        for medium, rich in MEDIA.items():
            with model:
                set_medium(model, rich, strict=strict)
                out[(name, medium)] = predict_essential(model)
            print(f"  {name} / {medium}: {int(out[(name, medium)]['predicted_essential'].sum())} "
                  "predicted essential", flush=True)
    return out


def grade_all(preds: dict, table, votes, genes=None) -> tuple:
    """Score every model x medium x experiment, optionally restricted to a fixed set of genes.

    `genes` exists because the models do not contain the same genes: iJR904 has 904, iML1515 has 1515.
    Scoring each model over its own gene list means a comparison between two models also changes which
    genes are being graded, so an apparent improvement could be a change of subject. Passing the genes
    common to every model holds that fixed.
    """
    rows, split_rows = [], []
    for (name, medium), pred in preds.items():
        call = pred["predicted_essential"]
        if genes is not None:
            call = call.reindex(genes).dropna().astype(bool)
        for src in SOURCES:
            rows.append({"model": name, "medium": medium, "experiment": src,
                         **score(call, table[src])})
            for subset, mask in (("agreed", votes["unanimous"]), ("contested", votes["contested"])):
                keep = table.index[mask.reindex(table.index).fillna(False)]
                split_rows.append({"model": name, "medium": medium, "experiment": src,
                                   "genes": subset, **score(call, table.loc[keep, src])})
    return pd.DataFrame(rows), pd.DataFrame(split_rows)


def main():
    RESULTS.mkdir(exist_ok=True)
    table = load()

    # ---------------------------------------------------------------- the experiments against each other
    agreement(table, SOURCES).to_csv(RESULTS / "agreement.csv", index=False)
    votes = consensus(table, SOURCES)
    votes.to_csv(RESULTS / "consensus.csv")
    called = votes[votes["n_essential"] > 0]
    print(f"genes called essential by at least one experiment: {len(called)}; "
          f"unanimous {int(called['unanimous'].sum())} ({called['unanimous'].mean():.0%}), "
          f"contested {int(called['contested'].sum())} ({called['contested'].mean():.0%})\n")

    # ---------------------------------------------------------------- grade every model
    preds = predictions()
    grades, split = grade_all(preds, table, votes)
    grades.to_csv(RESULTS / "grades.csv", index=False)
    split.to_csv(RESULTS / "split.csv", index=False)

    print(grades[grades["model"] == BEFORE]
          .pivot(index="medium", columns="experiment", values="mcc").round(3).to_string(), "\n")

    # ------------------------------------------------- comparing models without also changing the genes
    # The four generations do not contain the same genes, so a straight comparison between them is
    # partly a comparison of gene lists. Grading them all on the genes they share separates the two.
    generations = [n for n in MODEL_FILES if n != AFTER]
    common = set.intersection(*(set(preds[(n, "rich")].index) for n in generations))
    common = sorted(common)
    common_grades, _ = grade_all({k: v for k, v in preds.items() if k[0] in generations},
                                 table, votes, genes=common)
    common_grades["genes"] = "common"
    common_grades.to_csv(RESULTS / "common_genes.csv", index=False)
    sizes = ", ".join("{}: {}".format(n, len(preds[(n, "rich")])) for n in generations)
    print(f"the four generations share {len(common)} genes ({sizes} in total). "
          "Graded on those shared genes only:")
    print(common_grades.pivot_table(index="medium", columns="model", values="mcc")
          .reindex(columns=generations).round(3).to_string())
    own = grades[grades["model"].isin(generations)].pivot_table(index="medium", columns="model",
                                                                values="mcc")
    print("\nand for comparison, each graded on its own genes:")
    print(own.reindex(columns=generations).round(3).to_string(), "\n")

    # ------------------------------------------------- how much of that is the medium each model ships
    # set_medium keeps uptakes a model's author left open. That is right for grading a published model
    # and wrong for attributing a difference between two of them, so the strict run drops them.
    strict_preds = predictions(names=tuple(generations), strict=True)
    strict_grades, _ = grade_all(strict_preds, table, votes, genes=common)
    strict_grades["genes"] = "common"
    strict_grades.to_csv(RESULTS / "strict_uptakes.csv", index=False)
    print("same shared genes, but every model on exactly the same uptakes and nothing else:")
    print(strict_grades.pivot_table(index="medium", columns="model", values="mcc")
          .reindex(columns=generations).round(3).to_string(), "\n")

    # ---------------------------------------------------------------- the grade depends on the grader
    print("same model, same medium, graded on genes the experiments agree about vs contest:")
    s = split[(split["model"] == BEFORE) & (split["medium"] == "rich")]
    print(s.pivot(index="experiment", columns="genes", values="mcc").round(3).to_string())
    c = s[s["genes"] == "contested"]
    print("  on contested genes the experiments do not even agree which way the model is wrong:")
    for _, r in c.iterrows():
        print(f"    {r['experiment']:7s} calls {r['observed_essential']:3.0f} of {r['n']:3.0f} essential; "
              f"model recall {r['recall']:.3f}, precision {r['precision']:.3f}")
    print()

    # ---------------------------------------------------------------- do the published corrections hold up?
    # The corrections in iML1515_corrected were derived against RB-TnSeq carbon-source fitness data
    # (Bernstein et al. 2023). None of the three essentiality experiments here was used to make them, so
    # this is an out-of-sample test by construction, with no need to trust anybody's discipline.
    corrections = []
    for medium in MEDIA:
        for src in SOURCES:
            before = grades.query("model == @BEFORE and medium == @medium and experiment == @src")
            after = grades.query("model == @AFTER and medium == @medium and experiment == @src")
            corrections.append({
                "medium": medium, "experiment": src,
                "mcc_before": float(before["mcc"].iloc[0]), "mcc_after": float(after["mcc"].iloc[0]),
                "delta_mcc": float(after["mcc"].iloc[0] - before["mcc"].iloc[0]),
                "errors_before": int(before["fp"].iloc[0] + before["fn"].iloc[0]),
                "errors_after": int(after["fp"].iloc[0] + after["fn"].iloc[0]),
            })
    corr = pd.DataFrame(corrections)
    corr.to_csv(RESULTS / "corrections.csv", index=False)
    print("published corrections, graded on experiments they were not fitted to:")
    print(corr.round(3).to_string(index=False), "\n")

    # the error count and the MCC move in opposite directions under minimal medium; this is why
    rec = reconcile(grades[grades["model"].isin([BEFORE, AFTER])], BEFORE, AFTER)
    rec.to_csv(RESULTS / "reconcile.csv", index=False)
    print("the two measures disagree because they price the two classes differently:")
    print(rec.round(3).to_string(index=False), "\n")

    # what the corrections actually are, so the transfer result can be read next to its cause
    diff = corrections_diff(load_model(BEFORE), load_model(AFTER))
    diff["gene_rules_changed"].to_csv(RESULTS / "corrections_diff.csv", index=False)
    n_changed = len(diff["gene_rules_changed"])
    print(f"the corrections consist of: {len(diff['reactions_added'])} added reactions "
          f"{diff['reactions_added']}, {len(diff['reactions_removed'])} removed, "
          f"{n_changed + diff['gene_rules_equivalent']} gene rules rewritten of which "
          f"{diff['gene_rules_equivalent']} are logically equivalent and {n_changed} are not, "
          f"{diff['bounds_changed']} changed bounds")
    print(diff["gene_rules_changed"][
        ["reaction", "change", "dropped", "newly_essential", "no_longer_essential"]
    ].to_string(index=False), "\n")

    # which genes changed, and whether the change is bigger than noise
    chg = changed_calls(preds[(BEFORE, "rich")]["predicted_essential"],
                        preds[(AFTER, "rich")]["predicted_essential"], table, SOURCES)
    chg.to_csv(RESULTS / "changed.csv")
    print(f"the corrections change {len(chg)} calls under rich medium:")
    print(chg[["gene", "before", "after", *SOURCES, "moved"]].to_string(), "\n")
    print(direction(chg, SOURCES).to_string(index=False), "\n")

    boot = []
    for medium in MEDIA:
        for src in SOURCES:
            boot.append({"medium": medium, "experiment": src,
                         **bootstrap_delta(preds[(BEFORE, medium)]["predicted_essential"],
                                           preds[(AFTER, medium)]["predicted_essential"], table[src])})
    bootstrap = pd.DataFrame(boot)
    bootstrap.to_csv(RESULTS / "bootstrap.csv", index=False)
    print("resampling genes, 90% interval on the change in MCC:")
    print(bootstrap.round(4).to_string(index=False), "\n")

    # ---------------------------------------------------------------- where the model's errors sit
    pred = preds[(BEFORE, "rich")]["predicted_essential"]
    joined = pd.DataFrame({"predicted": pred}).join(votes[["n_sources", "contested", "unanimous"]])
    joined = joined.join(table[list(SOURCES)])
    wrong = pd.Series(False, index=joined.index)
    for src in SOURCES:
        wrong |= (joined["predicted"] != joined[src]) & joined[src].notna()
    share = joined.loc[wrong & joined["n_sources"].gt(1), "contested"].mean()
    print(f"of the genes iML1515 gets wrong against at least one experiment, "
          f"{share:.0%} are genes the experiments themselves disagree about")

    # ---------------------------------------------------------------- threshold sweep
    rows = []
    for threshold in (0.001, 0.01, 0.05, 0.1, 0.2):
        for medium in MEDIA:
            p = preds[(BEFORE, medium)]
            call = p["growth_fraction"] < threshold
            for src in SOURCES:
                rows.append({"threshold": threshold, "medium": medium, "experiment": src,
                             **score(call, table[src])})
    pd.DataFrame(rows).to_csv(RESULTS / "threshold.csv", index=False)
    print("\nwritten to results/")


if __name__ == "__main__":
    main()
