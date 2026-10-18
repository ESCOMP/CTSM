"""
Unit test for ctsm.machine_learning.export_torchscript
"""

import unittest
from ctsm.machine_learning.export_torchscript import SampleIdentityModel, ConstantModel, export_model_to_torchscript, HAS_TORCH

if HAS_TORCH:
    import torch


class TestExportTorchscript(unittest.TestCase):
    def test_sample_identity_model(self):
        """Test instantiation and identity behavior of SampleIdentityModel."""
        model = SampleIdentityModel()
        self.assertIsNotNone(model)
        if HAS_TORCH:
            x = torch.tensor([1.0, -2.5, 3.14])
            out = model(x)
            self.assertTrue(torch.allclose(out, x))

    def test_constant_model(self):
        """Test instantiation and constant behavior of ConstantModel."""
        model = ConstantModel(constant_value=2.0)
        self.assertIsNotNone(model)
        self.assertEqual(model.constant_value, 2.0)
        if HAS_TORCH:
            x = torch.tensor([0.5, 1.5, -3.0])
            out = model(x)
            expected = torch.full_like(x, 2.0)
            self.assertTrue(torch.allclose(out, expected))


if __name__ == "__main__":
    unittest.main()
