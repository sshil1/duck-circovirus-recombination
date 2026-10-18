#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate cavphylo
cd ~/dcv
for b in blockL blockR; do
  FastTree -nt -gtr -gamma rec3_$b.fasta > rec3_$b.nwk 2> rec3_$b.fasttree.log
  echo "$b tree done: $(tail -1 rec3_$b.fasttree.log)"
done
python3 -I step21_cmp.py
echo FINISHED; read
