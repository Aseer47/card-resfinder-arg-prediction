"""
arg_resistance_tools
=====================

Standalone, reusable tools extracted from the ARG classifier project
(Ruthbah, Sadi, Jahan & Adib, 2026 — "Sequence Composition-Based Prediction
of Antibiotic Resistance Genes"). Two independent capabilities are provided:

1. Resistance mechanism classification + per-class accuracy breakdown
   (with Wilson score 95% confidence intervals), so that classifier
   performance can be reported as a function of resistance *mechanism*
   (beta-lactamase, efflux pump, etc.) rather than as a single pooled
   number.

2. Cross-database sequence redundancy auditing via CD-HIT, so that a
   "cross-database generalization" claim can be checked against how much
   the two databases actually overlap at the protein sequence level
   before any performance number is interpreted.

Install
-------
    pip install git+https://github.com/Aseer47/card-resfinder-arg-prediction

Requires a working `cd-hit` binary on PATH for the redundancy-audit
functions (e.g. `conda install -c bioconda cd-hit` or
`apt-get install cd-hit`). The mechanism-classification functions have
no external dependencies beyond BioPython and NumPy/pandas.

Quick usage
-----------
    from arg_resistance_tools import classify_mechanism, mechanism_accuracy_breakdown

    mechanism_labels = [classify_mechanism(h) for h in card_headers]
    result = mechanism_accuracy_breakdown(y_true, y_pred, mechanism_labels)
    print(result.summary())
    result.plot()

    from arg_resistance_tools import cross_database_redundancy

    redundancy = cross_database_redundancy(
        "card_positive.fasta", "resfinder_translated.fasta",
        thresholds=(90, 70, 50),
    )
    print(redundancy.summary())
    redundancy.plot()
"""

from __future__ import annotations

import math
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from Bio import SeqIO
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "arg_resistance_tools requires BioPython. Install with: pip install biopython"
    ) from exc


# ---------------------------------------------------------------------------
# 1. Resistance mechanism classification
# ---------------------------------------------------------------------------

#: Default keyword -> mechanism-class mapping, as used in the ARG paper.
#: Matching is case-insensitive and checks whether any keyword appears as a
#: substring of the gene symbol / FASTA header.
DEFAULT_MECHANISM_KEYWORDS: Dict[str, Tuple[str, ...]] = {
    "Beta-lactamase": ("BLA", "AMPC", "OXA", "KPC", "NDM", "CTX", "SHV", "TEM"),
    "Aminoglycoside-modifying enzyme": ("AAC", "APH", "ANT"),
    "Tetracycline resistance": ("TET",),
    "Efflux pump component": ("MEX", "ACR", "OQX", "ADE", "MDS", "MTR", "SME"),
    "RNA polymerase target mutation": ("ROB", "RPO"),
}

#: Label used when no keyword matches.
UNCLASSIFIED_LABEL = "Other/unclassified ARG"


def classify_mechanism(
    header_or_gene_symbol: str,
    keyword_map: Optional[Dict[str, Tuple[str, ...]]] = None,
    unclassified_label: str = UNCLASSIFIED_LABEL,
) -> str:
    """Classify a CARD gene name or FASTA header into a resistance mechanism class.

    Parameters
    ----------
    header_or_gene_symbol:
        A CARD FASTA header (e.g. ``gb|ACCESSION|ARO:3002999|blaTEM-1``) or a
        bare gene symbol (e.g. ``blaTEM-1``). Matching is substring-based and
        case-insensitive, so either form works.
    keyword_map:
        Optional override of the default mechanism-class keyword mapping, in
        case a different gene nomenclature or an expanded set of mechanism
        classes is needed. Defaults to ``DEFAULT_MECHANISM_KEYWORDS``.
    unclassified_label:
        Label returned when no keyword in ``keyword_map`` matches.

    Returns
    -------
    str
        The matched mechanism class name, or ``unclassified_label``.

    Notes
    -----
    Keyword matching is first-match-wins in the order the keyword_map is
    defined. A gene symbol matching multiple keyword lists (rare, but
    possible with short substrings) is assigned to whichever class is
    checked first.
    """
    keyword_map = keyword_map or DEFAULT_MECHANISM_KEYWORDS
    text = header_or_gene_symbol.upper()
    for mechanism_class, keywords in keyword_map.items():
        if any(keyword in text for keyword in keywords):
            return mechanism_class
    return unclassified_label


def wilson_score_interval(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """Compute a Wilson score 95% confidence interval for a binomial proportion.

    Preferred over the standard Wald interval for small samples or
    proportions near 0% or 100%, where the normal approximation underlying
    the Wald interval breaks down (e.g. a reported 100% accuracy on n=1).

    Parameters
    ----------
    k:
        Number of successes (e.g. correctly classified sequences).
    n:
        Number of trials (e.g. total sequences in that category).
    z:
        Z-score for the desired confidence level. Defaults to 1.96 (95% CI).

    Returns
    -------
    (lower, upper):
        Interval bounds as percentages (0-100), not proportions.
    """
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half_width = (z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / denom
    lo = max(0.0, center - half_width)
    hi = min(1.0, center + half_width)
    return (lo * 100, hi * 100)


@dataclass
class MechanismAccuracyResult:
    """Per-mechanism-class accuracy breakdown with Wilson score confidence intervals."""

    table: pd.DataFrame  # columns: mechanism_class, n, accuracy_pct, ci_low, ci_high
    worst_class: str = field(init=False)
    best_class: str = field(init=False)

    def __post_init__(self) -> None:
        sorted_table = self.table.sort_values("accuracy_pct")
        self.worst_class = sorted_table.iloc[0]["mechanism_class"]
        self.best_class = sorted_table.iloc[-1]["mechanism_class"]

    def summary(self) -> str:
        lines = ["Resistance mechanism class accuracy breakdown:", ""]
        for _, row in self.table.sort_values("n", ascending=False).iterrows():
            lines.append(
                f"  {row['mechanism_class']:<35s} "
                f"n={int(row['n']):>4d}  "
                f"acc={row['accuracy_pct']:>5.1f}%  "
                f"95% CI=[{row['ci_low']:.1f}%, {row['ci_high']:.1f}%]"
            )
        lines.append("")
        lines.append(
            f"Best-performing class:  {self.best_class} "
            f"({self.table.set_index('mechanism_class').loc[self.best_class, 'accuracy_pct']:.1f}%)"
        )
        lines.append(
            f"Worst-performing class: {self.worst_class} "
            f"({self.table.set_index('mechanism_class').loc[self.worst_class, 'accuracy_pct']:.1f}%)"
        )
        return "\n".join(lines)

    def plot(self, ax=None):
        """Bar chart of per-class accuracy with Wilson CI error bars.

        Requires matplotlib (not a hard dependency of this module).
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 5))

        sorted_table = self.table.sort_values("accuracy_pct")
        y_pos = np.arange(len(sorted_table))
        acc = sorted_table["accuracy_pct"].values
        lo_err = acc - sorted_table["ci_low"].values
        hi_err = sorted_table["ci_high"].values - acc

        ax.barh(y_pos, acc, xerr=[lo_err, hi_err], capsize=4, color="#4C72B0")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(
            [
                f"{row['mechanism_class']} (n={int(row['n'])})"
                for _, row in sorted_table.iterrows()
            ]
        )
        ax.set_xlabel("Accuracy (%)")
        ax.set_xlim(0, 105)
        ax.set_title("Classification accuracy by resistance mechanism class")
        plt.tight_layout()
        return ax


def mechanism_accuracy_breakdown(
    y_true: Sequence[int],
    y_pred: Sequence[int],
    mechanism_labels: Sequence[str],
    positive_only: bool = True,
) -> MechanismAccuracyResult:
    """Compute per-mechanism-class classification accuracy with 95% CIs.

    Parameters
    ----------
    y_true:
        Ground-truth binary labels (1 = ARG, 0 = non-ARG), aligned with
        ``y_pred`` and ``mechanism_labels``.
    y_pred:
        Predicted binary labels, same length and order as ``y_true``.
    mechanism_labels:
        Mechanism class label for each sequence (e.g. from
        ``classify_mechanism``), same length and order as ``y_true``.
    positive_only:
        If True (default, matching the ARG paper's convention), restricts
        the breakdown to sequences where ``y_true == 1``, since mechanism
        class is only meaningful for true resistance genes. If False, all
        sequences are included regardless of true label.

    Returns
    -------
    MechanismAccuracyResult
        Object with a `.table` DataFrame, `.summary()` string method, and
        `.plot()` visualization method.
    """
    df = pd.DataFrame(
        {"y_true": y_true, "y_pred": y_pred, "mechanism_class": mechanism_labels}
    )
    if positive_only:
        df = df[df["y_true"] == 1]

    rows = []
    for mechanism_class, group in df.groupby("mechanism_class"):
        n = len(group)
        k = int((group["y_pred"] == group["y_true"]).sum())
        acc_pct = 100.0 * k / n if n else float("nan")
        ci_low, ci_high = wilson_score_interval(k, n)
        rows.append(
            {
                "mechanism_class": mechanism_class,
                "n": n,
                "accuracy_pct": acc_pct,
                "ci_low": ci_low,
                "ci_high": ci_high,
            }
        )

    table = pd.DataFrame(rows).reset_index(drop=True)
    return MechanismAccuracyResult(table=table)


# ---------------------------------------------------------------------------
# 2. Cross-database redundancy audit (CD-HIT)
# ---------------------------------------------------------------------------

def _parse_cdhit_clstr(clstr_path: str, prefix_a: str, prefix_b: str) -> List[List[str]]:
    """Parse a CD-HIT .clstr file into a list of clusters (lists of sequence IDs).

    Uses plain substring search rather than a regex on the sequence ID, since
    CD-HIT IDs can contain characters (parentheses, apostrophes — common in
    gene symbols such as ``aac(6')-Ib``) that break naive ID-extraction
    regexes. This is the fixed parser from the ARG paper's debugging process.
    """
    clusters: List[List[str]] = []
    current: List[str] = []

    with open(clstr_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = []
            elif prefix_a in line or prefix_b in line:
                for prefix in (prefix_a, prefix_b):
                    if prefix in line:
                        start = line.index(f">{prefix}") + 1
                        end = line.index("...", start)
                        current.append(line[start:end])
    if current:
        clusters.append(current)
    return clusters


@dataclass
class RedundancyResult:
    """Cross-database redundancy audit results across one or more CD-HIT identity thresholds."""

    table: pd.DataFrame  # columns: threshold, n_total_b, n_redundant, pct_redundant, n_nonredundant
    nonredundant_ids: Dict[int, List[str]]  # threshold -> list of non-redundant seq IDs from database B

    def summary(self) -> str:
        lines = ["Cross-database redundancy audit (CD-HIT):", ""]
        for _, row in self.table.iterrows():
            lines.append(
                f"  {int(row['threshold'])}% identity: "
                f"{int(row['n_redundant'])}/{int(row['n_total_b'])} redundant "
                f"({row['pct_redundant']:.1f}%), "
                f"{int(row['n_nonredundant'])} non-redundant remaining"
            )
        min_redundancy = self.table["pct_redundant"].min()
        lines.append("")
        if min_redundancy > 90:
            lines.append(
                f"WARNING: even at the lowest tested threshold, {min_redundancy:.1f}% "
                "of database B is redundant with database A. Cross-database "
                "generalization claims between these two databases should be "
                "interpreted with this redundancy in mind."
            )
        return "\n".join(lines)

    def plot(self, ax=None):
        """Bar chart of % redundancy across tested identity thresholds."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(6, 4))

        ax.bar(
            [str(int(t)) for t in self.table["threshold"]],
            self.table["pct_redundant"],
            color="#C44E52",
        )
        ax.set_xlabel("CD-HIT identity threshold (%)")
        ax.set_ylabel("Redundant with database A (%)")
        ax.set_ylim(0, 100)
        ax.set_title("Cross-database sequence redundancy")
        for i, pct in enumerate(self.table["pct_redundant"]):
            ax.text(i, pct + 1.5, f"{pct:.1f}%", ha="center")
        plt.tight_layout()
        return ax


def cross_database_redundancy(
    fasta_a: str,
    fasta_b: str,
    thresholds: Iterable[int] = (90, 70, 50),
    cdhit_path: str = "cd-hit",
    memory_mb: int = 2000,
    label_a: str = "A",
    label_b: str = "B",
) -> RedundancyResult:
    """Audit sequence-level redundancy between two protein FASTA databases via CD-HIT.

    For each identity threshold, combines both FASTA files, runs CD-HIT, and
    reports what fraction of database B's sequences fall into a cluster that
    also contains a database A sequence (i.e. are "redundant" with A). This
    should be run, and the result inspected, *before* interpreting any
    "trained on A, tested on B" cross-database performance claim.

    Parameters
    ----------
    fasta_a, fasta_b:
        Paths to the two protein FASTA files to compare (e.g. a training
        database and an independent test database). Sequences must already
        be in protein (not nucleotide) form; translate nucleotide sequences
        before calling this function.
    thresholds:
        CD-HIT sequence identity thresholds to test, as percentages.
        Word length is chosen automatically per CD-HIT's own guidance
        (5 for thresholds >= 60%, 4 for 60% > t >= 50%, 3 below that).
    cdhit_path:
        Path to the `cd-hit` binary, if not on PATH.
    memory_mb:
        Memory limit passed to CD-HIT's `-M` flag (0 = unlimited).
    label_a, label_b:
        Short labels used as FASTA ID prefixes to distinguish the two
        databases' sequences after combining them. Must not collide with
        any existing ID content.

    Returns
    -------
    RedundancyResult
        Object with a `.table` DataFrame, `.summary()` string method, a
        `.plot()` visualization method, and `.nonredundant_ids` giving the
        actual non-redundant database-B sequence IDs at each threshold
        (useful for building a genuinely non-redundant test set).
    """
    records_a = list(SeqIO.parse(fasta_a, "fasta"))
    records_b = list(SeqIO.parse(fasta_b, "fasta"))
    prefix_a, prefix_b = f"{label_a}__", f"{label_b}__"

    rows = []
    nonredundant_ids: Dict[int, List[str]] = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        combined_path = os.path.join(tmpdir, "combined.fasta")
        with open(combined_path, "w") as out:
            for rec in records_a:
                out.write(f">{prefix_a}{rec.id}\n{str(rec.seq)}\n")
            for rec in records_b:
                out.write(f">{prefix_b}{rec.id}\n{str(rec.seq)}\n")

        for threshold in thresholds:
            if threshold >= 60:
                word_length = 5
            elif threshold >= 50:
                word_length = 3
            else:
                raise ValueError(
                    f"CD-HIT identity threshold {threshold}% is below 50%, "
                    "which CD-HIT does not support directly."
                )

            out_path = os.path.join(tmpdir, f"out_{threshold}")
            subprocess.run(
                [
                    cdhit_path,
                    "-i", combined_path,
                    "-o", out_path,
                    "-c", str(threshold / 100),
                    "-n", str(word_length),
                    "-M", str(memory_mb),
                    "-d", "0",
                ],
                check=True,
                capture_output=True,
            )

            clusters = _parse_cdhit_clstr(f"{out_path}.clstr", prefix_a, prefix_b)

            redundant_b_ids = set()
            for cluster in clusters:
                has_a = any(seq_id.startswith(prefix_a) for seq_id in cluster)
                b_ids_in_cluster = [
                    seq_id[len(prefix_b):]
                    for seq_id in cluster
                    if seq_id.startswith(prefix_b)
                ]
                if has_a and b_ids_in_cluster:
                    redundant_b_ids.update(b_ids_in_cluster)

            all_b_ids = {rec.id for rec in records_b}
            nonredundant = sorted(all_b_ids - redundant_b_ids)
            nonredundant_ids[threshold] = nonredundant

            n_total_b = len(records_b)
            n_redundant = len(redundant_b_ids)
            rows.append(
                {
                    "threshold": threshold,
                    "n_total_b": n_total_b,
                    "n_redundant": n_redundant,
                    "pct_redundant": 100.0 * n_redundant / n_total_b if n_total_b else float("nan"),
                    "n_nonredundant": len(nonredundant),
                }
            )

    table = pd.DataFrame(rows)
    return RedundancyResult(table=table, nonredundant_ids=nonredundant_ids)


__all__ = [
    "DEFAULT_MECHANISM_KEYWORDS",
    "UNCLASSIFIED_LABEL",
    "classify_mechanism",
    "wilson_score_interval",
    "MechanismAccuracyResult",
    "mechanism_accuracy_breakdown",
    "RedundancyResult",
    "cross_database_redundancy",
]
