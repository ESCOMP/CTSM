#!/bin/bash
set -euo pipefail

# Initialize conda for this non-interactive shell
eval "$(/opt/conda/bin/conda shell.bash hook)"

# Accept conda Terms of Service for default channels
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create the ctsm_pylib conda environment
cd /workspace
./py_env_create -y -o

# Auto-activate ctsm_pylib in new terminal sessions
echo 'conda activate ctsm_pylib' >> ~/.zshrc
echo 'conda activate ctsm_pylib' >> ~/.bashrc

# Verify
conda activate ctsm_pylib
python --version
pylint --version
black --version
echo "CTSM conda environment setup complete."
