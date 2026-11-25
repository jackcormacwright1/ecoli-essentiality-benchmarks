"""The charts, from the files written by 02_grade.py.

Usage:  PYTHONPATH=src python scripts/03_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"

INK, MUTED, GRID = "#1d1d1b", "#6b6a66", "#e1e0d9"
BLUE, ORANGE, GREEN, PURPLE, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#b4b2ab"
# Experiments keep one colour across every figure, assigned in a fixed order and never recycled.
COLOUR = {"keio": BLUE, "tradis": ORANGE, "pec": GREEN}
NAME = {"keio": "Keio deletions", "tradis": "TraDIS transposon", "pec": "PEC curated"}
ORDER = ["keio", "tradis", "pec"]
MODELS = ["iJR904", "iAF1260", "iJO1366", "iML1515"]
# The pathway each changed gene belongs to. Named rather than guessed from the prefix, because the
# grouping is the mechanism: the two added exchanges feed folate and NAD, and those are what move.
PATHWAY = {
    "folE": "folate", "folB": "folate", "folK": "folate", "folP": "folate", "folC": "folate",
    "pabA": "folate", "pabB": "folate",
    "nadD": "NAD", "nadE": "NAD",
    "aroE": "aromatic", "aroK": "aromatic",
    "metE": "methionine", "thrA": "threonine",
}

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 160, "lines.linewidth": 2.0,
})


def save(fig, name):
    fig.savefig(RESULTS / name, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def figure_disagreement(consensus, agreement):
    """How much the experiments disagree, and where the disagreement sits."""
    # Restricted to the genes on which all three have an opinion, so the three bars are one population
    # and sum to the total quoted beside them.
    full = consensus[consensus["n_sources"] == 3]
    counts = [int((full["n_essential"] == k).sum()) for k in (3, 2, 1)]
    total = sum(counts)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.2),
                             gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.42})

    bars = axes[0].bar(range(3), counts, color=[GREEN, ORANGE, PURPLE], width=0.62)
    axes[0].bar_label(bars, fmt="%d", padding=3, fontsize=9.5, color=INK)
    axes[0].set(xticks=range(3), ylabel="genes", title="Called essential by…")
    axes[0].set_xticklabels(["all three\nexperiments", "exactly two", "exactly one"], fontsize=9)
    axes[0].set_ylim(0, max(counts) * 1.22)
    axes[0].annotate(f"{total - counts[0]} of these {total} genes\nare contested by somebody",
                     (0.97, 0.88), xycoords="axes fraction", ha="right", va="top", fontsize=9,
                     color=MUTED)

    # The disagreement is not symmetric: each pair splits by which source calls the gene essential, and
    # each bar is labelled with that source's name, so no legend is needed to read identity.
    pairs = agreement[agreement["source_a"].isin(ORDER) & agreement["source_b"].isin(ORDER)]
    y = np.arange(len(pairs))
    widest = max(pairs["essential_in_a_only"].max(), pairs["essential_in_b_only"].max())
    for yi, (_, r) in enumerate(pairs.iterrows()):
        for offset, src, n in ((-0.19, r["source_a"], r["essential_in_a_only"]),
                               (0.19, r["source_b"], r["essential_in_b_only"])):
            axes[1].barh(yi + offset, n, height=0.34, color=COLOUR[src])
            axes[1].text(n + widest * 0.02, yi + offset, f"{src} only: {int(n)}", va="center",
                         fontsize=8.5, color=MUTED)
    axes[1].set_yticks(y, [f"{a}\nvs {b}" for a, b in zip(pairs["source_a"], pairs["source_b"])],
                       fontsize=8.5)
    axes[1].set(xlabel="genes one source calls essential and the other does not",
                xlim=(0, widest * 1.42), title="Agreement is 96–99% overall — but only because\n"
                                               "19 genes in 20 are non-essential and easy")
    axes[1].grid(axis="y", visible=False)
    fig.suptitle("Three experiments, one organism, and no agreed answer to grade against",
                 x=0.01, ha="left", y=1.04)
    save(fig, "disagreement.png")


def figure_medium(grades):
    """The medium moves the grade more than fifteen years of model development does."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    for ax, medium, title in zip(axes, ("minimal", "rich"),
                                 ("M9 glucose (the model's default)", "Rich medium (the experiment's)")):
        g = grades[(grades["medium"] == medium) & (grades["model"].isin(MODELS))]
        for src in ORDER:
            s = g[g["experiment"] == src].set_index("model").reindex(MODELS)
            ax.plot(range(len(MODELS)), s["mcc"], color=COLOUR[src], marker="o", ms=6, label=NAME[src])
        ax.set(xticks=range(len(MODELS)), title=title, ylim=(0, 0.85))
        ax.set_xticklabels(MODELS, fontsize=8.5)
    axes[0].set_ylabel("Matthews correlation")
    # Two of the three lines sit within 0.01 of each other at the right edge, so end-of-line names would
    # overprint. The legend carries identity instead.
    axes[0].legend(fontsize=8.5, loc="lower right")
    fig.suptitle("Getting the medium right is worth about as much as fifteen years of model development",
                 x=0.01, ha="left", y=1.03)
    save(fig, "medium.png")


def figure_grader(split):
    """The same model scores differently depending on which experiment grades it, and on which genes."""
    s = split[(split["model"] == "iML1515") & (split["medium"] == "rich")]
    agreed = s[s["genes"] == "agreed"].set_index("experiment").reindex(ORDER)
    contested = s[s["genes"] == "contested"].set_index("experiment").reindex(ORDER)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"wspace": 0.38})
    x = np.arange(len(ORDER))
    b1 = axes[0].bar(x - 0.19, agreed["mcc"], width=0.34, color=[COLOUR[s] for s in ORDER])
    b2 = axes[0].bar(x + 0.19, contested["mcc"], width=0.34,
                     color=[COLOUR[s] for s in ORDER], alpha=0.42)
    axes[0].bar_label(b1, fmt="%.2f", padding=3, fontsize=9, color=INK)
    axes[0].bar_label(b2, fmt="%.2f", padding=3, fontsize=9, color=MUTED)
    axes[0].axhline(0, color=INK, lw=1.0)
    axes[0].set(xticks=x, ylabel="Matthews correlation", ylim=(-0.12, 0.92),
                title="One model, two kinds of gene")
    axes[0].set_xticklabels([NAME[s].replace(" ", "\n", 1) for s in ORDER], fontsize=8.5)
    # The hue already carries the experiment, so the key explains the shading only, in neutral grey.
    shading = [matplotlib.patches.Patch(facecolor=MUTED, label="genes the experiments agree about"),
               matplotlib.patches.Patch(facecolor=MUTED, alpha=0.42, label="genes they contest")]
    axes[0].legend(handles=shading, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.14))

    # on the contested genes the three experiments disagree in opposite directions, so the model's
    # apparent failure mode is a property of the grader rather than of the model
    share = contested["observed_essential"] / contested["n"]
    axes[1].barh(x, share, color=[COLOUR[s] for s in ORDER], height=0.55)
    for xi, (sh, n_ess, n) in enumerate(zip(share, contested["observed_essential"], contested["n"])):
        axes[1].text(sh + 0.02, xi, f"{int(n_ess)} of {int(n)}", va="center", fontsize=9, color=MUTED)
    axes[1].set_yticks(x, [NAME[s] for s in ORDER], fontsize=9)
    axes[1].set(xlim=(0, 1.18), xlabel="share of the contested genes this experiment calls essential",
                title="And they contest them in opposite directions")
    axes[1].grid(axis="y", visible=False)
    fig.suptitle("The model's accuracy is partly a property of whoever is grading it",
                 x=0.01, ha="left", y=1.02)
    save(fig, "grader.png")


def figure_transfer(bootstrap, changed):
    """The published corrections, graded on the three experiments they were not fitted to."""
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5),
                             gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.46})

    b = bootstrap.copy()
    b["label"] = [f"{m}, {NAME[e].split()[0]}" for m, e in zip(b["medium"], b["experiment"])]
    b = b.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(b))
    axes[0].hlines(y, b["delta_mcc_lo"], b["delta_mcc_hi"],
                   color=[COLOUR[e] for e in b["experiment"]], lw=3.0, alpha=0.45)
    axes[0].plot(b["delta_mcc_mean"], y, "o", ms=8,
                 color="none", markerfacecolor="none")
    for yi, (m, e) in enumerate(zip(b["delta_mcc_mean"], b["experiment"])):
        axes[0].plot([m], [yi], "o", ms=7.5, color=COLOUR[e])
    axes[0].axvline(0, color=INK, lw=1.2)
    axes[0].set_yticks(y, b["label"], fontsize=8.5)
    axes[0].set(xlabel="change in Matthews correlation (negative = worse)\n"
                       "point and 90% interval over resampled genes",
                title="Six out of six point the wrong way,\nthree of six clear the interval")
    axes[0].grid(axis="y", visible=False)

    # the thirteen genes whose call changed, which is where the mechanism shows
    c = changed.copy()
    c["pathway"] = [PATHWAY.get(g, "other") for g in c["gene"]]
    c = c.sort_values(["pathway", "gene"])
    y = np.arange(len(c))
    tone = {"toward": GREEN, "away": ORANGE, "unclear": GREY}
    axes[1].barh(y, [1 if m == "toward" else -1 for m in c["moved"]],
                 color=[tone[m] for m in c["moved"]], height=0.6)
    axes[1].set_yticks(y, [f"{g}  ({p})" for g, p in zip(c["gene"], c["pathway"])], fontsize=8.5)
    axes[1].axvline(0, color=INK, lw=1.2)
    axes[1].set(xlim=(-1.6, 1.6), xticks=[-1, 1], title="The thirteen calls that changed")
    axes[1].set_xticklabels(["away from\nthe experiments", "toward\nthe experiments"], fontsize=8.5)
    axes[1].grid(axis="y", visible=False)
    fig.suptitle("Corrections fitted to carbon-source fitness data do not transfer to essentiality",
                 x=0.01, ha="left", y=1.12)
    save(fig, "transfer.png")


def main():
    consensus = pd.read_csv(RESULTS / "consensus.csv", index_col=0)
    agreement = pd.read_csv(RESULTS / "agreement.csv")
    grades = pd.read_csv(RESULTS / "grades.csv")
    split = pd.read_csv(RESULTS / "split.csv")
    bootstrap = pd.read_csv(RESULTS / "bootstrap.csv")
    changed = pd.read_csv(RESULTS / "changed.csv")

    figure_disagreement(consensus, agreement)
    figure_medium(grades)
    figure_grader(split)
    figure_transfer(bootstrap, changed)


if __name__ == "__main__":
    main()
