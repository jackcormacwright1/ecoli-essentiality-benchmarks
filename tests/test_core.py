"""Tests.

Three kinds, in increasing order of what they are worth:

  1. arithmetic — `score`, `agreement`, `consensus`, `reconcile` on hand-made input where the answer can
     be worked out by eye. These catch the ordinary mistakes.
  2. behaviour — `set_medium` and `predict_essential` against cobrapy's own small textbook model, which
     ships with the package, so this runs without the 50MB of SBML the real analysis needs.
  3. the claims — the numbers quoted in the README, re-derived from the written results. These are the
     ones that matter: they fail if a later change moves a headline figure without anybody noticing.
     They skip rather than fail when the results have not been generated, so CI stays green on a clone
     that does not carry the models.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from eeb.essentiality import agreement, consensus, load
from eeb.model import M9_GLUCOSE, RICH_SUPPLEMENTS, predict_essential, score, set_medium
from eeb.transfer import bootstrap_delta, changed_calls, corrections_diff, direction, reconcile

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
MODELS = ROOT / "data" / "external" / "models"
needs_results = pytest.mark.skipif(not (RESULTS / "grades.csv").exists(),
                                   reason="run scripts/02_grade.py first")
needs_models = pytest.mark.skipif(not (MODELS / "iML1515.xml").exists(),
                                  reason="run scripts/01_fetch.py --download first")


@pytest.fixture(scope="module")
def textbook():
    import cobra
    return cobra.io.load_model("textbook")


# ----------------------------------------------------------------------------------- arithmetic
def test_score_counts_a_known_confusion_matrix():
    predicted = pd.Series([True, True, False, False, True], index=list("abcde"))
    observed = pd.Series([True, False, False, True, True], index=list("abcde"))
    s = score(predicted, observed)
    assert (s["tp"], s["fp"], s["tn"], s["fn"]) == (2, 1, 1, 1)
    assert s["precision"] == pytest.approx(2 / 3)
    assert s["recall"] == pytest.approx(2 / 3)
    assert s["accuracy"] == pytest.approx(0.6)


def test_score_ignores_genes_where_either_side_is_silent():
    predicted = pd.Series([True, False, True], index=list("abc"))
    observed = pd.Series([True, pd.NA, False], index=list("abc"), dtype="boolean")
    assert score(predicted, observed)["n"] == 2


def test_accuracy_flatters_and_mcc_does_not():
    """Calling nothing essential scores well on accuracy and zero on MCC. This asymmetry is the whole
    reason MCC is the headline measure, so it is worth pinning."""
    observed = pd.Series([True] * 5 + [False] * 95)
    nothing = pd.Series([False] * 100)
    s = score(nothing, observed)
    assert s["accuracy"] == pytest.approx(0.95)
    assert pd.isna(s["mcc"]) or s["mcc"] == 0


def test_agreement_counts_each_direction_separately():
    table = pd.DataFrame({
        "a": pd.array([True, True, False, False], dtype="boolean"),
        "b": pd.array([True, False, True, False], dtype="boolean"),
    })
    row = agreement(table, ("a", "b")).iloc[0]
    assert row["n_compared"] == 4 and row["agree"] == 2
    assert row["essential_in_a_only"] == 1 and row["essential_in_b_only"] == 1


def test_consensus_needs_more_than_one_opinion_to_call_a_gene_agreed():
    table = pd.DataFrame({
        "gene": ["only_one", "all_three", "split"],
        "a": pd.array([True, True, True], dtype="boolean"),
        "b": pd.array([None, True, False], dtype="boolean"),
        "c": pd.array([None, True, False], dtype="boolean"),
    })
    out = consensus(table, ("a", "b", "c"))
    assert not out.loc[0, "unanimous"] and not out.loc[0, "contested"]  # one voice is not a consensus
    assert out.loc[1, "unanimous"] and out.loc[1, "consensus_essential"]
    assert out.loc[2, "contested"]


def test_reconcile_shows_the_trade_behind_a_falling_error_count():
    """A change that sheds more false positives than it loses true positives reduces the error count.
    `reconcile` has to report both sides of that trade, because the error count alone hides it."""
    grades = pd.DataFrame([
        {"model": "before", "medium": "m", "experiment": "x", "n": 100, "observed_essential": 10,
         "tp": 8, "fp": 20, "fn": 2, "recall": 0.8, "precision": 0.286, "mcc": 0.40},
        {"model": "after", "medium": "m", "experiment": "x", "n": 100, "observed_essential": 10,
         "tp": 5, "fp": 10, "fn": 5, "recall": 0.5, "precision": 0.333, "mcc": 0.35},
    ])
    r = reconcile(grades, "before", "after").iloc[0]
    assert r["tp_lost"] == 3 and r["fp_shed"] == 10
    assert r["delta_errors"] == -7  # errors fall
    assert r["delta_mcc"] < 0       # and the model is still worse


# ----------------------------------------------------------------------------------- behaviour
def test_set_medium_closes_everything_it_did_not_open(textbook):
    with textbook as model:
        applied = set_medium(model, rich=False)
        open_uptakes = {r.id for r in model.exchanges if r.lower_bound < 0}
    assert open_uptakes == set(applied)
    assert all(ex in M9_GLUCOSE or ex in RICH_SUPPLEMENTS for ex in applied)


def test_rich_medium_opens_more_than_minimal(textbook):
    with textbook as model:
        minimal = set_medium(model, rich=False)
    with textbook as model:
        rich = set_medium(model, rich=True)
    assert set(minimal) < set(rich)


def test_set_medium_keeps_an_uptake_the_author_deliberately_left_open(textbook):
    """iML1515_corrected arrives with almost every uptake shut and three left open on purpose. Those
    three are part of the correction, so putting every model on one base must not discard them."""
    with textbook as model:
        model.reactions.EX_ac_e.lower_bound = -5.0  # stand-in for an author's extra
        applied = set_medium(model, rich=False)
    assert applied.get("EX_ac_e") == -5.0


def test_predict_essential_refuses_a_model_that_cannot_grow(textbook):
    """The published corrected model is infeasible as distributed. Silently scoring it returns zero
    essential genes and a plausible-looking accuracy, which is the worst possible failure."""
    with textbook as model:
        for exchange in model.exchanges:
            exchange.lower_bound = 0.0
        with pytest.raises(ValueError, match="does not grow"):
            predict_essential(model)


def test_predict_essential_finds_some_essential_genes(textbook):
    with textbook as model:
        set_medium(model, rich=False)
        out = predict_essential(model)
    assert 0 < out["predicted_essential"].sum() < len(out)
    assert out["growth_fraction"].between(0, 1.001).all()


def test_changed_calls_reports_direction_against_the_majority():
    before = pd.Series([True, False, True], index=list("abc"))
    after = pd.Series([False, True, True], index=list("abc"))
    observed = pd.DataFrame({
        "gene": ["g1", "g2", "g3"],
        "keio": pd.array([False, False, True], dtype="boolean"),
        "pec": pd.array([False, False, True], dtype="boolean"),
    }, index=list("abc"))
    chg = changed_calls(before, after, observed, ("keio", "pec"))
    assert list(chg.index) == ["a", "b"]          # c did not change
    assert chg.loc["a", "moved"] == "toward"      # True -> False, and both say non-essential
    assert chg.loc["b", "moved"] == "away"
    d = direction(chg, ("keio", "pec")).set_index("experiment")
    assert d.loc["keio", "now_right"] == 1 and d.loc["keio", "now_wrong"] == 1


def test_corrections_diff_ignores_rebracketing_but_not_real_changes():
    """Most of the edits in the published corrected model only re-bracket an unchanged rule. Counting
    those as changes would treble the apparent size of the correction.

    Built from stubs rather than cobra models because cobra normalises rule text on assignment, which
    would quietly erase the very distinction under test.
    """
    class R:
        def __init__(self, rid, rule="", bounds=(0.0, 1000.0)):
            self.id, self.gene_reaction_rule, self.bounds = rid, rule, bounds

    class M:
        def __init__(self, reactions):
            self.reactions = reactions

    before = M([R("KEPT", "b0001 and b0002"), R("REBRACKET", "(b0003 and b0004) or b0005"),
                R("TIGHTENED", "b0006 or b0007"), R("REBOUNDED", "b0008", (0.0, 1000.0))])
    after = M([R("KEPT", "b0001 and b0002"), R("REBRACKET", "b0003 and b0004 or b0005"),
               R("TIGHTENED", "b0006"), R("REBOUNDED", "b0008", (-7.0, 1000.0)),
               R("EX_added_c", "", (-1000.0, 1000.0))])

    d = corrections_diff(before, after)
    assert d["reactions_added"] == ["EX_added_c"] and d["reactions_removed"] == []
    assert d["gene_rules_rebracketed"] == 1
    assert list(d["gene_rules_changed"]["reaction"]) == ["TIGHTENED"]
    assert d["gene_rules_changed"].iloc[0]["dropped"] == "b0007"
    assert d["gene_rules_changed"].iloc[0]["added"] == ""
    assert d["bounds_changed"] == 1


@needs_models
def test_the_published_corrections_are_what_the_readme_says_they_are():
    """The README explains the transfer failure by what the corrections contain. If the model file ever
    changes under us, that explanation stops being true and this fails rather than going quietly stale."""
    from eeb.model import load_model
    d = corrections_diff(load_model("iML1515"), load_model("iML1515_corrected"))
    assert d["reactions_added"] == ["EX_nad_c", "EX_thf_c"]
    assert d["reactions_removed"] == []
    rules = d["gene_rules_changed"].set_index("reaction")
    assert len(rules) == 9
    # every substantive rule change drops an isozyme and adds nothing
    assert (rules["added"] == "").all()
    assert rules.loc["METS", "rule_after"] == "b3829"     # metE left alone when metH is dropped
    assert rules.loc["SHKK", "rule_after"] == "b3390"     # aroK left alone when aroL is dropped


def test_bootstrap_delta_is_zero_when_nothing_changed():
    pred = pd.Series([True, False] * 25, index=range(50))
    observed = pd.Series([True, False] * 25, index=range(50), dtype="boolean")
    out = bootstrap_delta(pred, pred, observed, draws=25)
    assert out["delta_mcc_mean"] == pytest.approx(0.0)
    assert out["n_genes"] == 50


# ----------------------------------------------------------------------------------- the claims
def test_the_experiments_disagree_about_half_the_genes_they_call_essential():
    table = load()
    votes = consensus(table, ("keio", "tradis", "pec"))
    called = votes[votes["n_essential"] > 0]
    assert len(called) == 402
    assert int(called["contested"].sum()) == 187


def test_tn5_never_calls_anything_essential():
    """The reason tn5 is loaded but kept out of the headline comparison. If this ever fails, the
    exclusion needs revisiting rather than leaving in place."""
    assert not load()["tn5"].fillna(False).any()


@needs_results
def test_the_corrections_do_not_transfer():
    corr = pd.read_csv(RESULTS / "corrections.csv")
    assert len(corr) == 6
    assert (corr["delta_mcc"] < 0).all(), "every out-of-sample test should still go the wrong way"


@needs_results
def test_the_error_count_and_the_mcc_disagree_under_minimal_medium():
    """The finding that made the write-up honest: raw errors fall while the model gets worse."""
    rec = pd.read_csv(RESULTS / "reconcile.csv").set_index(["medium", "experiment"])
    assert rec.loc[("minimal", "keio"), "delta_errors"] < 0
    assert rec.loc[("minimal", "keio"), "delta_mcc"] < 0
    assert (rec["tp_lost"] > 0).all()


@needs_results
def test_the_grade_depends_on_which_genes_it_is_taken_over():
    split = pd.read_csv(RESULTS / "split.csv")
    s = split[(split["model"] == "iML1515") & (split["medium"] == "rich")]
    agreed = s[s["genes"] == "agreed"]["mcc"]
    contested = s[s["genes"] == "contested"]["mcc"]
    assert agreed.min() > 0.75           # on agreed genes the three graders are interchangeable
    assert agreed.max() - agreed.min() < 0.01
    assert contested.max() < 0.45        # on contested genes they are not
