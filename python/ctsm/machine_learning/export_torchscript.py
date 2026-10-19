#!/usr/bin/env python3
"""
Tool and library to create toy models and export PyTorch models to TorchScript.

As a tool, this creates toy models for testing. As a library, it converts a
PyTorch nn.Module to a TorchScript file.
"""

import argparse
import logging
import sys

try:
    import torch  # pylint: disable=import-error
    from torch import nn  # pylint: disable=import-error

    HAS_TORCH = True
except ImportError:
    # Allows ctsm python test discovery (run_ctsm_py_tests) to succeed
    # on environments without PyTorch installed.
    torch = None
    nn = None
    HAS_TORCH = False

logger = logging.getLogger(__name__)


if HAS_TORCH:

    class SampleIdentityModel(nn.Module):  # pylint: disable=too-few-public-methods
        """Sample Identity Model used for testing and demonstration.

        Returns input unchanged.
        """

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            """Forward pass returning the input tensor unchanged."""
            return x * 1.0

    class ConstantModel(nn.Module):  # pylint: disable=too-few-public-methods
        """Sample Constant Model used for testing and demonstration.

        Returns tensor of ones matching input shape, multiplied by a constant.
        """

        def __init__(self, constant_value: float = 1.0):
            super().__init__()
            self.constant_value = constant_value

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            """Forward pass returning a constant tensor of the input shape."""
            return torch.ones_like(x) * self.constant_value

else:

    class SampleIdentityModel:  # pylint: disable=too-few-public-methods
        """Fallback SampleIdentityModel placeholder when PyTorch is unavailable."""

    class ConstantModel:  # pylint: disable=too-few-public-methods
        """Fallback ConstantModel placeholder when PyTorch is unavailable."""

        def __init__(self, constant_value: float = 1.0):
            self.constant_value = constant_value


def export_model_to_torchscript(model, output_file, method="script", example_input=None):
    """Export a PyTorch nn.Module to a TorchScript file.

    :param model: PyTorch nn.Module instance
    :param output_file: Path to output .pt file
    :param method: 'script' or 'trace'
    :param example_input: Tensor input required if method is 'trace'
    """
    if not HAS_TORCH:
        raise RuntimeError(
            "PyTorch is not installed in the active environment. Please activate ctsm_pylib."
        )

    model.eval()
    if method == "script":
        scripted_model = torch.jit.script(model)
    elif method == "trace":
        if example_input is None:
            raise ValueError("example_input is required when method='trace'")
        scripted_model = torch.jit.trace(model, example_input)
    else:
        raise ValueError(f"Unknown export method: {method}. Must be 'script' or 'trace'.")

    scripted_model.save(output_file)
    logger.info("Successfully exported TorchScript model to: %s", output_file)
    return output_file


def main():
    """Command-line interface for exporting sample models to TorchScript."""
    parser = argparse.ArgumentParser(
        description="Export PyTorch model to TorchScript .pt format for FTorch."
    )
    parser.add_argument(
        "-o",
        "--output",
        default="constant_model.pt",
        help="Output .pt file path (default: constant_model.pt)",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=["script", "trace"],
        default="script",
        help="Export method (default: script)",
    )
    parser.add_argument(
        "--model-type",
        choices=["constant", "identity"],
        default="constant",
        help="Type of model to export (default: constant)",
    )
    parser.add_argument(
        "--constant-val",
        type=float,
        default=1.0,
        help="Value for constant model (default: 1.0)",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if not HAS_TORCH:
        print(
            "ERROR: PyTorch is not available. Ensure PyTorch is installed in your "
            "python environment (ctsm_pylib).",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.model_type == "constant":
        model = ConstantModel(constant_value=args.constant_val)
    else:
        model = SampleIdentityModel()

    export_model_to_torchscript(model, args.output, method=args.method)
    print(f"Exported {args.model_type} TorchScript model to: {args.output}")


if __name__ == "__main__":
    main()
