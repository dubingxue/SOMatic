"""Feature construction for SOMatic."""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdPartialCharges


def charge_context(smiles: str, n_iter: int = 12) -> np.ndarray:
    """Return six Gasteiger descriptors plus six missingness flags per heavy atom.

    Descriptors are atomic charge, attached-H charge, one-hop sum and mean, and
    exactly-two-hop sum and mean. Missing continuous values are emitted as NaN;
    :class:`SOMaticModel` standardizes them on training data, imputes zero, and
    appends the flags returned here.
    """
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.AddHs(base)
    rdPartialCharges.ComputeGasteigerCharges(mol, nIter=n_iter)
    heavy = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    result = []
    for atom in heavy:
        q = _charge(atom)
        hydrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
        hq = sum((_charge(h) for h in hydrogens), 0.0) if hydrogens else 0.0
        one = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        one_ids = {n.GetIdx() for n in one}
        two = {
            n2.GetIdx(): n2
            for n1 in one
            for n2 in n1.GetNeighbors()
            if n2.GetAtomicNum() > 1 and n2.GetIdx() != atom.GetIdx()
            and n2.GetIdx() not in one_ids
        }
        one_q = np.asarray([_charge(n) for n in one], dtype=float)
        two_q = np.asarray([_charge(n) for n in two.values()], dtype=float)
        values = np.asarray([
            q, hq, _sum_or_nan(one_q), _mean_or_nan(one_q),
            _sum_or_nan(two_q), _mean_or_nan(two_q),
        ])
        result.append(np.r_[values, ~np.isfinite(values)])
    return np.asarray(result, dtype=np.float32)


def _charge(atom: Chem.Atom) -> float:
    try:
        value = float(atom.GetProp("_GasteigerCharge"))
    except (KeyError, ValueError):
        return float("nan")
    return value if np.isfinite(value) else float("nan")


def _sum_or_nan(values: np.ndarray) -> float:
    return float(np.sum(values)) if len(values) and np.all(np.isfinite(values)) else float("nan")


def _mean_or_nan(values: np.ndarray) -> float:
    return float(np.mean(values)) if len(values) and np.all(np.isfinite(values)) else float("nan")
