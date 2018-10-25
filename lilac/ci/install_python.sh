#!/usr/bin/env bash

set -e
set -x

# Install miniconda
wget --quiet http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /usr/src/miniconda.sh
bash /usr/src/miniconda.sh -b -p /usr/local/miniconda
conda update conda --yes
conda clean -tipy
conda config --set always_yes yes --set changeps1 no
conda --version

conda install -c conda-forge cmake>=3
