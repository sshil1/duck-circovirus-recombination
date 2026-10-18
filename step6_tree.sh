#!/usr/bin/env bash
# STEP 6 - ML tree. RUN INSIDE TMUX (IQ-TREE with 1000 UFBoot on ~443 taxa takes tens of minutes).
# Support is quoted as SH-aLRT/UFBoot; 'well supported' = SH-aLRT>=80 AND UFBoot>=95.
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh; conda activate cavphylo
T=${THREADS:-16}
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true); [ -z "$IQ" ] && { echo "IQ-TREE not found"; exit 1; }
FT=$(command -v FastTree || command -v fasttree || true); [ -z "$FT" ] && { echo "FastTree not found"; exit 1; }
"$FT" -nt -gtr -gamma tree_aln.fasta > dcv_ft.nwk 2> dcv_ft.log && echo "FastTree done -> dcv_ft.nwk"
"$IQ" -s tree_aln.fasta -m MFP -B 1000 -alrt 1000 -T "$T" --prefix dcv_iqtree -redo
echo "IQ-TREE done -> dcv_iqtree.treefile  (model/lnL in dcv_iqtree.iqtree)"
