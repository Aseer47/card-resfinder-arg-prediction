"""
example_usage.py
=================

Minimal, runnable example of both tools in arg_resistance_tools.py:

1. Resistance mechanism classification + per-class accuracy breakdown
   (synthetic toy data; no external binary required).
2. Cross-database redundancy audit via CD-HIT (synthetic toy FASTA files;
   requires a `cd-hit` binary on PATH — this part is skipped with a
   friendly message if CD-HIT is not installed).

Run with:
    python example_usage.py
"""

import shutil
import tempfile
from pathlib import Path

from arg_resistance_tools import (
    classify_mechanism,
    mechanism_accuracy_breakdown,
    cross_database_redundancy,
)


def example_mechanism_breakdown() -> None:
    print("=" * 70)
    print("Example 1: resistance mechanism classification + accuracy breakdown")
    print("=" * 70)

    # Toy CARD-style headers, a mix of mechanism classes.
    headers = (
        ["gb|AAA{:05d}|ARO:3000000|blaTEM-1".format(i) for i in range(20)]
        + ["gb|AAA{:05d}|ARO:3000001|aac(6')-Ib".format(i) for i in range(10)]
        + ["gb|AAA{:05d}|ARO:3000002|tetA".format(i) for i in range(8)]
        + ["gb|AAA{:05d}|ARO:3000003|mexF".format(i) for i in range(6)]
        + ["gb|AAA{:05d}|ARO:3000004|rpoB".format(i) for i in range(1)]
    )
    mechanism_labels = [classify_mechanism(h) for h in headers]

    # Toy ground truth (all true ARGs) and toy predictions: the classifier
    # does well on beta-lactamases and aminoglycosides, poorly on the
    # structurally distinct efflux pump class, mirroring the paper's result.
    y_true = [1] * len(headers)
    y_pred = (
        [1] * 20          # beta-lactamase: 20/20 correct
        + [1] * 9 + [0]   # aminoglycoside: 9/10 correct
        + [1] * 6 + [0] * 2  # tetracycline: 6/8 correct
        + [1] * 2 + [0] * 4  # efflux pump: 2/6 correct
        + [1]             # RNA polymerase: 1/1 correct
    )

    result = mechanism_accuracy_breakdown(y_true, y_pred, mechanism_labels)
    print(result.summary())
    print()


def example_redundancy_audit() -> None:
    print("=" * 70)
    print("Example 2: cross-database redundancy audit (requires cd-hit)")
    print("=" * 70)

    if shutil.which("cd-hit") is None:
        print(
            "cd-hit binary not found on PATH — skipping this example.\n"
            "Install with: conda install -c bioconda cd-hit  (or apt-get install cd-hit)\n"
        )
        return

    # Toy protein FASTA: database A has 3 sequences, database B has 4, two
    # of which are near-identical to sequences in A (simulated redundancy)
    # and two of which are genuinely distinct.
    seq_shared = "MKTLLLTLVVVTIVCLDLGYT" * 3  # arbitrary repeated motif, >50 aa
    seq_distinct_1 = "GATTACAGATTACAGATTACA" * 5
    seq_distinct_2 = "WVWVWVWVWVWVWVWVWVWVW" * 5

    fasta_a = ">A_seq1\n{0}\n>A_seq2\n{1}\n>A_seq3\n{2}\n".format(
        seq_shared, seq_distinct_1[:60], seq_distinct_2[:60]
    )
    fasta_b = (
        ">B_seq1_redundant\n{0}\n"
        ">B_seq2_redundant\n{0}\n"
        ">B_seq3_novel\n{1}\n"
        ">B_seq4_novel\n{2}\n"
    ).format(seq_shared, seq_distinct_1[5:65], seq_distinct_2[5:65])

    with tempfile.TemporaryDirectory() as tmpdir:
        path_a = Path(tmpdir) / "db_a.fasta"
        path_b = Path(tmpdir) / "db_b.fasta"
        path_a.write_text(fasta_a)
        path_b.write_text(fasta_b)

        redundancy = cross_database_redundancy(
            str(path_a), str(path_b), thresholds=(90,), label_a="A", label_b="B"
        )
        print(redundancy.summary())
        print()
        print("Non-redundant database-B IDs at 90% identity:",
              redundancy.nonredundant_ids[90])


if __name__ == "__main__":
    example_mechanism_breakdown()
    example_redundancy_audit()
