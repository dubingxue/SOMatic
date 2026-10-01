import unittest
import numpy as np

from somatic.features import charge_context
from somatic.model import SOMaticModel, entity_balanced_weights


class FrameworkTest(unittest.TestCase):
    def test_feature_shape_and_fit(self):
        q = charge_context("CCO")
        self.assertEqual(q.shape, (3, 12))
        rng = np.random.default_rng(13)
        x = rng.normal(size=(6, 768)).astype(np.float32)
        charges = np.vstack([q, q])
        y = np.array([0, 1, 0, 0, 0, 1])
        entities = ["a"] * 3 + ["b"] * 3
        records = ["ra"] * 3 + ["rb"] * 3
        weights = entity_balanced_weights(entities, records, [3] * 6)
        model = SOMaticModel(13).fit(x, charges, y, weights)
        scores = model.predict_proba(x, charges)
        self.assertEqual(scores.shape, (6,))
        self.assertTrue(np.all((scores >= 0) & (scores <= 1)))
        self.assertAlmostEqual(float(weights.sum()), 1.0)


if __name__ == "__main__":
    unittest.main()
