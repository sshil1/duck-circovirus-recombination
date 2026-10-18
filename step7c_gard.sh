#!/usr/bin/env bash
# STEP 7c - GARD: locate recombination breakpoints on the block-aware subsample. RUN INSIDE TMUX (can take hours).
# GTR + 3-class GDD rate variation; 'Faster' mode (approximate but standard for >~50 taxa).
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh; conda activate hyphy
CPU=${CPU:-16}
hyphy gard --alignment gard_input.fasta --type nucleotide --rv GDD --rate-classes 3 --mode Faster \
           --output gard_result.json --output-lf gard_fit.bf CPU=$CPU 2>&1 | tee gard.log
echo "GARD done -> gard_result.json"
