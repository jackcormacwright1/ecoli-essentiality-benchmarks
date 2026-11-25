"""Comparing two gene-protein-reaction rules by what they mean.

A GPR rule is a Boolean expression over gene identifiers using only `and` and `or` — never `not`. Such
an expression is *monotone*, and monotone Boolean functions have an exact canonical form: the set of
**minimal** gene sets that are sufficient to enable the reaction. `(A and B) or A` and `A` have the same
canonical form; `A and B` and `A or B` do not. Two rules mean the same thing if and only if their
canonical forms are equal.

This matters because it is the difference that drives an essentiality prediction. `A and B` makes both
genes essential — knock out either and the reaction stops. `A or B` makes neither essential, because each
covers for the other. A comparison that looks only at which genes *appear* in a rule cannot tell those
two apart, and would report a complex being downgraded to a pair of isozymes as no change at all.

The construction is a small recursive-descent parser into disjunctive normal form, followed by
absorption: a term is dropped when another term is a subset of it, since the smaller set already
suffices. The result is the set of minimal sufficient gene sets, which is unique.
"""

from __future__ import annotations

import re

TOKEN = re.compile(r"\s*(\(|\)|\band\b|\bor\b|[A-Za-z0-9_.\-]+)", re.IGNORECASE)

# A rule with many independent `or` branches joined by `and` can expand combinatorially in DNF. Real GPR
# rules do not approach this, but a guard keeps a pathological input from hanging the analysis.
MAX_TERMS = 20000


class GPRTooComplex(ValueError):
    """Raised when a rule's disjunctive normal form exceeds MAX_TERMS."""


def _tokenise(rule: str) -> list[str]:
    tokens, position = [], 0
    while position < len(rule):
        match = TOKEN.match(rule, position)
        if match is None:
            if rule[position:].strip() == "":
                break
            raise ValueError(f"cannot parse gene rule near {rule[position:position + 20]!r}")
        token = match.group(1)
        tokens.append(token.lower() if token.lower() in ("and", "or") else token)
        position = match.end()
    return tokens


def _absorb(terms: set[frozenset[str]]) -> frozenset[frozenset[str]]:
    """Drop any term that a smaller term already covers, leaving only the minimal sufficient sets."""
    minimal = []
    for term in sorted(terms, key=len):
        if not any(kept <= term for kept in minimal):
            minimal.append(term)
    return frozenset(minimal)


def minimal_gene_sets(rule: str) -> frozenset[frozenset[str]]:
    """The canonical form of a GPR rule: the minimal sets of genes sufficient to enable the reaction.

    An empty rule returns the empty set, meaning no gene enables the reaction (a spontaneous or orphan
    reaction), which is distinct from a rule satisfied by no genes at all.
    """
    tokens = _tokenise(rule)
    if not tokens:
        return frozenset()
    position = 0

    def parse_or() -> set[frozenset[str]]:
        nonlocal position
        terms = parse_and()
        while position < len(tokens) and tokens[position] == "or":
            position += 1
            terms = terms | parse_and()
            if len(terms) > MAX_TERMS:
                raise GPRTooComplex(rule)
        return terms

    def parse_and() -> set[frozenset[str]]:
        nonlocal position
        terms = parse_atom()
        while position < len(tokens) and tokens[position] == "and":
            position += 1
            right = parse_atom()
            terms = {a | b for a in terms for b in right}
            if len(terms) > MAX_TERMS:
                raise GPRTooComplex(rule)
        return terms

    def parse_atom() -> set[frozenset[str]]:
        nonlocal position
        if position >= len(tokens):
            raise ValueError(f"gene rule ends unexpectedly: {rule!r}")
        token = tokens[position]
        if token == "(":
            position += 1
            inner = parse_or()
            if position >= len(tokens) or tokens[position] != ")":
                raise ValueError(f"unbalanced parentheses in gene rule: {rule!r}")
            position += 1
            return inner
        if token in ("and", "or", ")"):
            raise ValueError(f"unexpected {token!r} in gene rule: {rule!r}")
        position += 1
        return {frozenset({token})}

    parsed = parse_or()
    if position != len(tokens):
        raise ValueError(f"trailing tokens in gene rule: {rule!r}")
    return _absorb(parsed)


def equivalent(a: str, b: str) -> bool:
    """Do two GPR rules describe the same Boolean function?"""
    return minimal_gene_sets(a) == minimal_gene_sets(b)


def essential_genes(rule: str) -> frozenset[str]:
    """The genes whose single deletion disables the reaction.

    A gene is essential to a rule when it appears in every minimal sufficient set: with it gone, no
    remaining combination suffices. This is what connects a rule change to an essentiality prediction,
    and it is the quantity a gene-set comparison silently discards.
    """
    sets = minimal_gene_sets(rule)
    if not sets:
        return frozenset()
    return frozenset.intersection(*sets)


def describe_change(before: str, after: str) -> str:
    """A short label for how a rule changed, in terms of what it does to essentiality."""
    if equivalent(before, after):
        return "equivalent"
    gained = essential_genes(after) - essential_genes(before)
    lost = essential_genes(before) - essential_genes(after)
    if gained and not lost:
        return "more essential"
    if lost and not gained:
        return "less essential"
    if gained and lost:
        return "both"
    return "restructured"
