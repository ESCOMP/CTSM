"""Unit test for ctsm.machine_learning.export_torchscript."""

import os
import tempfile
import unittest

from ctsm.machine_learning.export_torchscript import (
    HAS_TORCH,
    ConstantModel,
    SampleIdentityModel,
    export_model_to_torchscript,
)

try:
    import torch  # pylint: disable=import-error
except ImportError:
    torch = None


class TestExportTorchscript(unittest.TestCase):
    """Unit tests for export_torchscript module and models."""

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

    def test_export_model(self):
        """Test export_model_to_torchscript behavior with or without torch."""
        if not HAS_TORCH:
            with self.assertRaises(RuntimeError):
                export_model_to_torchscript(None, "dummy.pt")
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_file = os.path.join(tmpdir, "model.pt")
                model = ConstantModel(constant_value=1.0)
                export_model_to_torchscript(model, output_file, method="script")
                self.assertTrue(os.path.isfile(output_file))


if __name__ == "__main__":
    unittest.main()
