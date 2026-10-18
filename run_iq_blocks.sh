#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate cavphylo
cd ~/dcv
V=$(iqtree --version 2>&1 | head -2 | tr '\n' ' ')
echo "$V"
if echo "$V" | grep -qE 'version (2|3)'; then BB=-B; else BB=-bb; fi
for b in blockL blockR; do
  iqtree -s rec3_$b.fasta -m GTR+G $BB 1000 -pre iq_$b -redo > iq_$b.stdout 2>&1
  echo "$b: $(grep -m1 'Log-likelihood of the tree' iq_$b.log)"
done
ls -la iq_block*.treefile
echo FINISHED; read
