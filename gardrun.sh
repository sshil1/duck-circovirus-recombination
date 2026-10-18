#!/usr/bin/env bash
# usage: gardrun.sh OUTDIR FASTA RV NP
set -euo pipefail
OUT=$1; FA=$2; RV=$3; NP=$4
E=""; MPI=""
for CAND in /opt/conda/envs/hyphympi /opt/conda/envs/hyphy; do
  for M in mpirun mpiexec orterun; do
    if [ -x $CAND/bin/HYPHYMPI ] && [ -x $CAND/bin/$M ]; then E=$CAND; MPI=$CAND/bin/$M; break 2; fi
  done
done
[ -n "$MPI" ] || { echo "NO MPI LAUNCHER FOUND"; ls /opt/conda/envs/*/bin/ 2>/dev/null | grep -i -E "mpi|orte" ; exit 1; }
mkdir -p ~/dcv/$OUT && cd ~/dcv/$OUT && cp ~/dcv/$FA .
echo "out=$OUT fasta=$FA rv=$RV ranks=$NP launcher=$MPI"
$MPI -np $NP $E/bin/HYPHYMPI "ENV=TOLERATE_NUMERICAL_ERRORS=1;" gard --alignment $FA --type nucleotide --rv $RV --mode Faster --output gard_result.json --output-lf gard_fit.bf 2>&1 | tee gard.log
