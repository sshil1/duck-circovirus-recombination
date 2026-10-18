#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate seq3env
cd ~/dcv/seq3
PT=/opt/conda/envs/seq3env/share/3seq/PVT.3SEQ.2017.700
yes Y | 3seq -full rec3_subset.fasta -ptable $PT 2>&1 | tee 3seq_full.log
echo FINISHED
ls -lt | head -8
read
