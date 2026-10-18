#!/usr/bin/env bash
# STEP 7b - separate conda env for HyPhy (keeps cavphylo untouched). Run in tmux.
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh
conda create -y -n hyphy -c conda-forge -c bioconda hyphy
conda activate hyphy
which hyphy; hyphy --version 2>&1 | head -2
echo "---- GARD options on THIS install (paste me this if step7c errors) ----"
hyphy gard --help 2>&1 | head -60
