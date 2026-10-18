#!/usr/bin/env bash
# STEP 0 - confirm the conda env has what steps 1-4 need. Reuses env 'cavphylo' from the CAV kit.
source /opt/conda/etc/profile.d/conda.sh && conda activate cavphylo || { echo "env cavphylo missing: run CAV kit 00_env.sh first"; exit 1; }
python3 -c "import Bio, numpy; print('biopython', Bio.__version__, '| numpy', numpy.__version__)" || conda install -y -n cavphylo -c conda-forge biopython numpy
command -v mafft >/dev/null && mafft --version 2>&1 | head -1 || conda install -y -n cavphylo -c bioconda mafft
echo "ready"
