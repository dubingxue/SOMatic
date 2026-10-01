"""Command-line trainer for SOMatic."""

from __future__ import annotations

import argparse
from pathlib import Path
import joblib
import numpy as np
import pandas as pd

from .features import charge_context
from .model import SOMaticModel, entity_balanced_weights


def atom_set(value: str) -> set[int]:
    return {int(x) for x in str(value).split(";") if x != ""}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True)
    p.add_argument("--atom-features", required=True)
    p.add_argument("--family", required=True, choices=["CYP", "UGT", "FMO", "AOX1"])
    p.add_argument("--output", required=True)
    p.add_argument("--seed", type=int, default=13, choices=[13, 37, 71])
    args = p.parse_args()

    records = pd.read_csv(args.records, dtype=str)
    records = records.loc[records.enzyme_family.eq(args.family)].copy()
    cache = np.load(args.atom_features)
    lookup = {(str(e), int(a)): x for e, a, x in zip(cache["entity_id"], cache["atom_index"], cache["features"])}

    xs, qs, ys, entities, record_ids, counts = [], [], [], [], [], []
    charge_cache = {}
    for row in records.itertuples(index=False):
        candidates, positives = atom_set(row.candidate_atoms), atom_set(row.positive_atoms)
        if not positives <= candidates:
            raise ValueError(f"positive atom outside candidate set: {row.record_id}")
        charge_cache.setdefault(row.entity_id, charge_context(row.smiles))
        for atom in sorted(candidates):
            xs.append(lookup[(row.entity_id, atom)])
            qs.append(charge_cache[row.entity_id][atom])
            ys.append(int(atom in positives))
            entities.append(row.entity_id)
            record_ids.append(row.record_id)
            counts.append(len(candidates))
    weights = entity_balanced_weights(entities, record_ids, counts)
    model = SOMaticModel(args.seed).fit(xs, qs, ys, weights)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "family": args.family, "seed": args.seed}, output)


if __name__ == "__main__":
    main()
