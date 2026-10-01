# SOMatic

SOMatic is an atom-level site-of-metabolism ranking framework. It combines a frozen 768-dimensional Uni-Mol2-84M atom representation with six local Gasteiger charge descriptors and six missing-value indicators, then fits a weighted ExtraTrees classifier.

This repository intentionally contains only the model framework and a traceable description of the data sources. Manuscript files, figures, private curated records, trained checkpoints, benchmark predictions, and external baseline code are not included.

## Model

- frozen Uni-Mol2-84M atom representation: 768 dimensions;
- charge context: atomic charge, attached-hydrogen charge, and sum/mean charge at one and exactly two graph hops;
- missing-value indicators: one per charge descriptor;
- predictor: 128-tree ExtraTrees, minimum leaf size 2, square-root feature sampling, Gini criterion, no bootstrap;
- technical seeds: 13, 37, and 71;
- fitting weight: equal total weight per molecular entity, divided equally over its records and candidate atoms;
- separate model per enzyme family; CYP subtype is not an input.

## Install

```bash
conda env create -f environment.yml
conda activate somatic
```

Uni-Mol2 weights are not redistributed here. Obtain them through the official `unimol_tools` package and its documented model-download mechanism.

## Input format

Prepare:

1. `records.csv`, following `examples/records.example.csv`;
2. `atom_features.npz`, containing arrays `features`, `entity_id`, and `atom_index`. `features` must contain one 768-dimensional frozen Uni-Mol2 vector per atom.

Candidate and positive atom indices are zero-based and semicolon separated. Every positive atom must be present in the candidate set.

## Train

```bash
python -m somatic.train \
  --records records.csv \
  --atom-features atom_features.npz \
  --family CYP \
  --output outputs/cyp_seed13.joblib \
  --seed 13
```

Repeat with seeds 37 and 71. Use `--family UGT`, `FMO`, or `AOX1` to fit the corresponding family-specific model.

## Test the framework

```bash
python -m unittest discover -s tests -v
```

The test uses synthetic frozen atom embeddings and does not download data or model weights.

## Data

See [DATA_SOURCES.md](DATA_SOURCES.md) for source publications, identifiers, access routes, and redistribution boundaries. Reproducing published numerical results additionally requires the corrected, versioned record table and split manifest described there; those source-derived data are not silently replaced with synthetic or reformatted labels.

## Scope

Unreported candidate atoms are proxy-zero training targets, not experimentally verified negatives. The public framework reproduces the declared feature construction, weighting, and estimator configuration; it does not claim that unavailable historical molecular coordinates can be regenerated exactly.
