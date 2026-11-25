# Methodology

What was decided, why, and what this design cannot tell you.

## The question

Metabolic model papers report gene essentiality accuracy against a single experiment, usually the Keio
collection. That number is treated as a property of the model. This repository asks how much of it is
instead a property of the experiment, the medium, and the growth threshold — three choices made before
the model is consulted at all.

## What a prediction is

A gene is predicted essential when, with the reactions its gene-protein-reaction rule permits switched
off, flux balance analysis can no longer carry flux through the biomass reaction above a threshold
fraction of wild-type growth. This is `cobra.flux_analysis.single_gene_deletion` with a threshold applied
afterwards, so the sweep over thresholds costs nothing extra.

Three choices sit inside that sentence.

### The medium

The largest one. A histidine biosynthesis gene is essential in minimal medium and dispensable in rich
medium, because in rich medium the cell imports histidine instead. iML1515's default is M9 with glucose.
The Keio collection was built on LB. Grading a minimal-medium model against rich-medium essentiality
scores every amino-acid biosynthesis gene as a false positive, and the resulting number says more about
the mismatch than about the model.

Both media are therefore run and both are reported. "Rich" is M9 glucose plus 38 exchange reactions —
twenty amino acids, six nucleosides, and the vitamins and cofactors a rich broth supplies — opened at
10 mmol gDW⁻¹ h⁻¹. It is a stand-in for LB, not a reconstruction of it: LB's actual composition varies
between labs and is not fully characterised. The gap between the two media is the finding; neither number
on its own should be quoted as *the* accuracy.

### Putting five models on one medium

The five models do not arrive on comparable media. `iML1515_model_adjusted_all_corrections.xml` is
**infeasible as distributed**: 24 of its uptake reactions are closed, glucose among them, and three are
left open (biotin, pantothenate, thiamine). The authors saved it part-way through an analysis that set
the medium itself.

This is a trap rather than an inconvenience. Run a deletion scan on it unchanged and the wild type does
not grow, every deletion also does not grow, no deletion reduces growth below a fraction of zero, and
every gene comes back dispensable — which looks like a result and scores a plausible-seeming accuracy
because 95% of genes really are dispensable.

`set_medium()` therefore puts every model on one base: a canonical M9 glucose taken from iML1515's own
default, **plus** whatever uptakes that model's author deliberately left open. The second half matters —
for the corrected model those three open uptakes are part of the correction, not an accident, and
discarding them would be grading a model nobody published. `predict_essential()` raises rather than
scoring a model that does not grow.

### Comparing models without also changing the genes

The four reconstructions contain 904, 1261, 1367 and 1516 genes. Scoring each over its own gene list
means that a comparison between two generations also changes which genes are being graded, so an
apparent improvement could be a change of subject rather than a change of quality. Every model is
therefore scored twice: once on its own genes, and once restricted to the 877 genes all four contain.

The two disagree, and the restricted version is the one that answers the question asked. On shared genes
the gain from iJR904 to iML1515 is concentrated almost entirely in the first step; the decade after 2007
contributes +0.008. The difference between the two rows also says something on its own — each
generation's additional genes are predicted less well than the shared core.

A second run puts every model on strictly identical uptakes (`set_medium(..., strict=True)`), discarding
any uptake a model's author left open, so that a difference between models cannot be carried by the
medium each one happens to ship. For the four published generations this changes nothing, because none
of them ships extras; the only model that does is the corrected one, where those three uptakes are part
of the correction. That is asserted by a test rather than assumed.

### The threshold

Deleting a gene rarely takes growth to exactly zero; it takes it to some fraction of wild type. Where the
cut is placed moves the essential set. It is set at 1% and swept over 0.1%–20% in `results/threshold.csv`.
Lethal deletions often return a small negative growth rate from the LP solver rather than exactly zero;
these are clipped to zero, which changes no call but keeps the reported fractions physical.

## The experiments

Five opinions are loaded; three are used in the headline comparison.

| Source | Basis | Used? |
|---|---|---|
| Keio (Baba et al. 2006) | no verified in-frame deletion mutant could be made | yes |
| TraDIS (Goodall et al. 2018) | gap in dense transposon insertion density | yes |
| PEC | curated call from the literature | yes |
| Tn5 | earlier MG1655 transposon screen, collated in the Keio supplement | no |
| genetic | classical genetic literature's call, collated in the Keio supplement | no |

Tn5 is excluded because it never calls a single gene essential in this table — it is a one-sided source
and cannot disagree with anything, so including it would inflate agreement statistics for free. A test
asserts this, so if the data ever changes the exclusion gets revisited rather than inherited. `genetic`
is excluded because it is a literature collation rather than one experiment, which makes it
non-independent of PEC in an unquantifiable way.

Calls are held as `True` / `False` / `NA`, with `NA` meaning the source declined to say. TraDIS is the
only source with an explicit "unclear" class; those genes become `NA` rather than being folded into
non-essential, which would otherwise silently convert an admission of uncertainty into a negative call.
Every comparison is restricted to genes where both sides have an opinion.

### Mapping

Everything is keyed on b-numbers. The Keio supplement carries them directly. TraDIS Table S4 is keyed by
gene name and is mapped through the Keio table's name→b-number correspondence; genes whose names do not
match are dropped rather than guessed at, which is why TraDIS comparisons have a smaller *n* (about 1,300
against about 1,500). The PEC export supplies a second opinion that agrees with the PEC column collated in
the Keio supplement wherever both speak, and fills the gaps where only one does.

### Why disagreement is expected

Technologies fail in different ways. A transposon screen miscalls a gene whose essential region is short,
because insertions elsewhere in it are tolerated, and miscalls the gene upstream of an essential one,
because the insertion knocks out a shared promoter. A deletion collection miscalls a gene whose mutant
survives only by acquiring a suppressor mutation, and inherits any error in library construction. So
where two sources disagree, neither is automatically the one that is wrong — and a model graded against
either is being graded partly on that source's failure modes.

This is the reason the analysis reports the split between genes the experiments agree on and genes they
contest, rather than picking a winner or taking a majority vote. A majority vote across three sources with
correlated failure modes manufactures a consensus that does not exist.

## The out-of-sample test

The corrections in `iML1515_corrected` were derived against RB-TnSeq carbon-source fitness data. None of
the three essentiality experiments used here contributed to them. This makes the test out-of-sample by
construction rather than by anybody's discipline about held-out data, which is the only form of the
claim worth making about somebody else's model.

Three things are reported together, because any one alone would be misleading:

**The summary statistic.** Matthews correlation before and after, for all six medium × experiment
combinations. MCC rather than accuracy because essential genes are 6–11% of the graded set and calling
everything dispensable scores about 93% accuracy.

**The reconciliation.** Under minimal medium the corrections reduce the raw error count while reducing
MCC. This is not a paradox and it is not a reason to prefer whichever measure flatters the conclusion:
counting errors prices a lost true positive and a lost false positive identically, and at 8% prevalence
they are not identical. `reconcile()` reports true positives lost and false positives shed side by side
so the trade is visible. A write-up that quoted only the error count could have called this result
positive.

**The bootstrap.** 600 resamples over genes, reported as a 90% interval. Three of the six tests exclude
zero on their own; all six point the same way. The interval is understated: genes within a pathway do not
vary independently, and the corrections act on pathways, so the effective sample size is smaller than the
gene count suggests. The six tests are also not independent of each other — one model, overlapping gene
sets, three experiments that agree about most genes. The defensible claim is a consistent small loss, not
a demonstrated large one.

**The diff.** What the corrections actually consist of, computed from the two model files rather than
described from the paper: two added intracellular exchange reactions, 37 rewritten gene rules of which 28
are logically equivalent to what they replaced and 9 are not, and 35 changed bounds.

Comparing gene rules needs the right test, and the obvious one is wrong. Asking whether the same genes
appear in a rule cannot distinguish `A and B` from `A or B` — identical gene sets, opposite meanings: a
complex in which both genes are essential, against isozymes in which neither is. Since a rule change
matters here only through its effect on essentiality, that is precisely the distinction that must survive
the comparison.

A GPR rule uses only `and` and `or`, never `not`, so it is a *monotone* Boolean function, and monotone
functions have an exact canonical form: the set of minimal gene sets sufficient to enable the reaction.
`src/eeb/gpr.py` parses each rule into disjunctive normal form and applies absorption — dropping any term
that a smaller term already covers — which yields that canonical form. Two rules mean the same thing if
and only if their canonical forms are equal. One real example needs the absorption step: `PFL` is written
as `(b0902 and b0903) or … or ((b0902 and b0903) and b2579)`, whose last term is redundant.

The same machinery gives the quantity that connects a rule to a prediction: a gene is essential to a rule
when it appears in **every** minimal sufficient set. Each of the nine real changes drops an isozyme and
makes exactly one gene newly essential.

## What this cannot tell you

**Which experiment is right.** Nothing here adjudicates a contested gene. The contested set is used as a
measurement of how much the grade depends on the grader, not as a problem to be solved.

**Whether the corrections are good corrections.** They were made for a different prediction task —
fitness across carbon sources — and this repository tests transfer to essentiality, which is not what
they were built for and not what they claimed. The precise result is that they reduced MCC on the three
essentiality benchmarks tested, under both media. That is evidence corrections are task-specific; it is
not evidence that the corrected model is worse at its own job, which is untested here. The README's
explanation of *why* they fail to transfer, via carry-over exchanges and isozyme removal, is read off the
model diff and is consistent with that reading.

**Whether rich medium is the right medium.** The rich medium here is a stand-in. A better analysis would
use a characterised LB composition, and the honest position is that the minimal–rich gap brackets the
answer rather than locating it.

**Anything about conditional essentiality.** Every gene gets one call per medium. Genes essential only
under specific stresses, or whose essentiality depends on growth rate, are outside this design entirely.

**Anything about non-metabolic genes.** A metabolic model has opinions only about the genes it contains,
roughly 1,500 of about 4,400. The large majority of essential genes in *E. coli* — ribosomal proteins,
DNA replication, cell division — are invisible to it. The accuracies here are over the metabolic subset
and are not comparable to whole-genome essentiality predictors.
