# arg-resistance-tools

Standalone utilities extracted from the ARG classifier project (Ruthbah, Sadi,
Jahan & Adib, 2026). Two independent capabilities, usable with any ARG
classifier or any pair of protein sequence databases, not just CARD/ResFinder.

## Install

```bash
pip install git+https://github.com/Aseer47/card-resfinder-arg-prediction
```

or, for the optional plotting methods:

```bash
pip install "arg-resistance-tools[plot] @ git+https://github.com/Aseer47/card-resfinder-arg-prediction"
```

## 1. Resistance mechanism classification + accuracy breakdown

Classify CARD gene names into resistance mechanism classes (beta-lactamase,
efflux pump, etc.) and compute per-class classification accuracy with Wilson
score 95% confidence intervals, rather than reporting a single pooled
accuracy number that can mask large differences between mechanism classes.

```python
from arg_resistance_tools import classify_mechanism, mechanism_accuracy_breakdown

mechanism_labels = [classify_mechanism(header) for header in card_headers]
result = mechanism_accuracy_breakdown(y_true, y_pred, mechanism_labels)

print(result.summary())
result.plot()  # requires matplotlib
```

## 2. Cross-database redundancy audit (CD-HIT)

Before interpreting a "trained on database A, tested on database B"
cross-database generalization claim, check how much the two databases
actually overlap at the protein sequence level.

```python
from arg_resistance_tools import cross_database_redundancy

redundancy = cross_database_redundancy(
    "card_positive.fasta", "resfinder_translated.fasta",
    thresholds=(90, 70, 50),
)

print(redundancy.summary())
redundancy.plot()  # requires matplotlib

# The genuinely non-redundant sequence IDs at each threshold, for building
# an actual held-out cross-database test set:
nonredundant_ids_at_90pct = redundancy.nonredundant_ids[90]
```

Requires a working `cd-hit` binary on PATH (e.g.
`conda install -c bioconda cd-hit` or `apt-get install cd-hit`).

## See also

See `example_usage.py` in this repository for a complete, runnable
end-to-end example of both tools.

## Citation

If you use this code, please cite:

Ruthbah, C.A., Sadi, T.H., Jahan, N.E.S., & Adib, A.N.M.T. (2026). Antibiotic
Resistance Gene Prediction from Sequence Composition: Performance,
Interpretability, and Database Redundancy. bioRxiv.

## License

MIT.
