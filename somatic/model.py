"""Weighted ExtraTrees implementation of the SOMatic predictor."""

from __future__ import annotations

import numpy as np
from sklearn.ensemble import ExtraTreesClassifier


def entity_balanced_weights(entity_ids, record_ids, candidate_counts) -> np.ndarray:
    """Equalize entity mass, then divide it over records and candidates."""
    entity_ids = np.asarray(entity_ids)
    record_ids = np.asarray(record_ids)
    candidate_counts = np.asarray(candidate_counts, dtype=float)
    n_entities = len(set(entity_ids.tolist()))
    # Each row represents one candidate. Count unique record IDs per entity.
    record_totals = {
        entity: len({record for e, record in zip(entity_ids, record_ids) if e == entity})
        for entity in set(entity_ids.tolist())
    }
    weights = np.asarray([
        1.0 / (n_entities * record_totals[e] * c)
        for e, c in zip(entity_ids, candidate_counts)
    ])
    return weights / weights.sum()


class SOMaticModel:
    """SOMatic charge preprocessing and ExtraTrees estimator."""

    def __init__(self, seed: int = 13):
        self.seed = seed
        self.charge_mean_: np.ndarray | None = None
        self.charge_scale_: np.ndarray | None = None
        self.estimator = ExtraTreesClassifier(
            n_estimators=128, min_samples_leaf=2, max_features="sqrt",
            criterion="gini", bootstrap=False, random_state=seed, n_jobs=-1,
        )

    def fit(self, atom_repr, charge_features, y, sample_weight):
        x = self._combine(atom_repr, charge_features, fit=True)
        self.estimator.fit(x, np.asarray(y, dtype=np.int8), sample_weight=sample_weight)
        return self

    def predict_proba(self, atom_repr, charge_features) -> np.ndarray:
        x = self._combine(atom_repr, charge_features, fit=False)
        return self.estimator.predict_proba(x)[:, 1]

    def _combine(self, atom_repr, charge_features, fit: bool) -> np.ndarray:
        atom_repr = np.asarray(atom_repr, dtype=np.float32)
        charge_features = np.asarray(charge_features, dtype=np.float32)
        if atom_repr.ndim != 2 or atom_repr.shape[1] != 768:
            raise ValueError("atom_repr must have shape (n_atoms, 768)")
        if charge_features.shape != (len(atom_repr), 12):
            raise ValueError("charge_features must have shape (n_atoms, 12)")
        continuous = charge_features[:, :6].copy()
        missing = charge_features[:, 6:].copy()
        if fit:
            self.charge_mean_ = np.nanmean(continuous, axis=0)
            self.charge_scale_ = np.nanstd(continuous, axis=0)
            self.charge_scale_[self.charge_scale_ == 0] = 1.0
        if self.charge_mean_ is None or self.charge_scale_ is None:
            raise RuntimeError("model has not been fitted")
        continuous = (continuous - self.charge_mean_) / self.charge_scale_
        continuous[~np.isfinite(continuous)] = 0.0
        return np.concatenate([atom_repr, continuous, missing], axis=1)
