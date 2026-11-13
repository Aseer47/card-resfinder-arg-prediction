CARD-ResFinder ARG Prediction

Interpretable ML framework for predicting antibiotic resistance genes (ARGs) from protein sequence composition, paired with a CD-HIT-based audit of how redundant two major ARG reference databases actually are before any "cross-database generalization" claim is made.

Key Findings

An XGBoost classifier trained on amino acid trigram frequencies and physicochemical descriptors separates ARGs from non-ARG proteins with near-ceiling accuracy (AUROC = 0.982, MCC = 0.886), confirmed as genuine signal rather than an artifact via a label-shuffling null baseline (AUROC = 0.508) and tight bootstrap 95% confidence intervals.

CARD and ResFinder, two of the field's most widely used ARG databases, are found to be 96.8-98.7% redundant at the protein sequence level, even at a CD-HIT identity threshold as permissive as 50%, indicating that cross-database generalization claims built on these two resources may largely be evaluating near-duplicate data.

On the small genuinely non-redundant residual of ResFinder sequences (40-100, depending on threshold), the classifier fails completely (0% classified as ARG at any threshold), and SHAP traces this failure to sequence length: non-redundant sequences average ~2.8x longer than the CARD sequences the model was trained on.

Classification accuracy varies sharply by resistance mechanism class: 100% for beta-lactamases (n=351, compact single-domain enzymes) versus 38.5% for efflux pump components (n=13, large multi-domain transporters), demonstrating that the same length/architecture axis explains both this in-distribution mechanism-class gap and the cross-database failure above.

Contents
├── ARG_prediction.ipynb        # Primary analysis pipeline (Colab-ready)
├── arg_resistance_tools.py     # Standalone reusable Python tool
├── example_usage.py            # Minimal usage example for the tool
├── setup.py                    # pip-installable package setup
├── README_tool.md              # Tool-specific documentation
├── results/                    # Figures, tables, robustness outputs
└── README.md
Data

Three public sequence databases, accessed directly via their official download endpoints. No new data was collected.

Class	Source	n	Role
Positive (ARG)	CARD (McArthur et al., 2013; Alcock et al., 2023)	4,835	Training + held-out test (80/20 split)
Negative (non-ARG)	UniProt/Swiss-Prot, BLAST-filtered against CARD	5,000	Training + held-out test (80/20 split)
Independent test	ResFinder (Bortolaia et al., 2020)	3,164 total; 40-100 non-redundant with CARD, depending on CD-HIT threshold	Cross-database generalization test
Reusable Tool

The resistance-mechanism classification and CD-HIT redundancy-audit methods are available as a standalone Python tool, installable directly from this repo:

bash
pip install git+https://github.com/Aseer47/card-resfinder-arg-prediction
python
from arg_resistance_tools import classify_mechanism, mechanism_accuracy_breakdown, cross_database_redundancy

mechanism_labels = [classify_mechanism(header) for header in card_headers]
result = mechanism_accuracy_breakdown(y_true, y_pred, mechanism_labels)
print(result.summary())
result.plot()

redundancy = cross_database_redundancy("card.fasta", "resfinder_translated.fasta", thresholds=(90, 70, 50))
print(redundancy.summary())
redundancy.plot()

See README_tool.md for full documentation.

Methods Summary
Amino acid trigram (3-mer) frequencies (8,000 features) plus five physicochemical descriptors (length, isoelectric point, aromaticity, instability index, GRAVY hydropathy)
Logistic Regression, Random Forest, and XGBoost evaluated on an 80/20 stratified train/test split
SHAP interpretability via TreeExplainer, with summary (beeswarm), categorical bar, and dependence plots
CD-HIT redundancy audit (90%/70%/50% identity) between CARD and ResFinder, performed before interpreting any cross-database result
Robustness checks: label-shuffling null baseline (10 permutations), bootstrap 95% confidence intervals (200 resamples)
Resistance mechanism class breakdown with Wilson score 95% confidence intervals
Running

Open ARG_prediction.ipynb in Google Colab and run all cells top to bottom. The CD-HIT/BLAST+ setup and BLAST filtering steps take several minutes on first run and are memory-constrained to fit the Colab free tier (-M 2000). All intermediate outputs checkpoint to Google Drive at every major stage, and final figures/tables save to results/.

Citation

If you use this code or the resistance-mechanism/redundancy-audit tool, please cite:

Ruthbah, C.A., Sadi, T.H., Jahan, N.E.S., & Adib, A.N.M.T. (2026). Sequence Composition-Based Prediction of Antibiotic Resistance Genes: Classifier Performance, Interpretability, and a Cross-Database Redundancy Audit. bioRxiv. (DOI to be added upon posting.)

Authors

Chowdhury Aseer Ruthbah, Talim Hossain Sadi, Nur E Shiratun Jahan, Abu Nayem Md. Tanzim Adib, BRAC University

License

MIT. CARD and UniProt sequences are redistributed under their original CC BY-SA 4.0 / CC BY 4.0 licenses.

Contact

Chowdhury Aseer Ruthbah: chowdhury.aseer.ruthbah@g.brac.ac.bd
