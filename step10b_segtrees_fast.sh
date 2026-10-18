#!/usr/bin/env bash
set -euo pipefail
set +u; source /opt/conda/etc/profile.d/conda.sh; conda activate cavphylo; set -u
cd ~/dcv/segB
ls seg?.fasta | sort -V | xargs -P 3 -I{} bash -c 'p=${1%.fasta}; echo "== $p $(date +%T)"; iqtree3 -s $1 -m GTR+F+G4 -B 1000 -bnni -alrt 1000 -nm 1000 -T 7 --prefix ${p}_f -redo -quiet; echo "done $p $(date +%T)"' _ {}
echo ALL_DONE $(date +%T)
