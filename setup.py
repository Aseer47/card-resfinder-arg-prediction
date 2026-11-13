from setuptools import setup

with open("README_tool.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="arg-resistance-tools",
    version="1.0.0",
    description=(
        "Resistance-mechanism classification with confidence intervals, and "
        "CD-HIT-based cross-database redundancy auditing, for antibiotic "
        "resistance gene (ARG) classifiers."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Chowdhury Aseer Ruthbah, Talim Hossain Sadi, Nur E Shiratun Jahan, Abu Nayem Md. Tanzim Adib",
    author_email="chowdhury.aseer.ruthbah@g.brac.ac.bd",
    url="https://github.com/Aseer47/card-resfinder-arg-prediction",
    license="MIT",
    py_modules=["arg_resistance_tools"],
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.21",
        "pandas>=1.3",
    ],
    extras_require={
        "plot": ["matplotlib>=3.4"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=(
        "antibiotic-resistance bioinformatics machine-learning shap "
        "cd-hit sequence-redundancy interpretability"
    ),
)
