#!/usr/bin/env bash
# STEP 5b - recombination screen (PHI + positive control). Seconds-to-minutes. Reads tree_aln.fasta.
# A negative PHI is reported ONLY together with the positive control (phi_validate.py): it bounds what the negative rules out.
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh; conda activate cavphylo
python3 -I phi_own.py tree_aln.fasta  | tee phi_test.txt
echo; echo "=== positive control ==="
python3 -I phi_validate.py            | tee phi_validation.txt
