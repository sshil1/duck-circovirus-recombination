#!/usr/bin/env bash
set -euo pipefail
for E in /opt/conda/envs/hyphympi /opt/conda/envs/hyphy; do
  [ -x $E/bin/HYPHYMPI ] && { ls $E/bin/mpirun >/dev/null 2>&1 || ls $E/bin/mpiexec >/dev/null 2>&1; } && break
done
MPI=$(ls $E/bin/mpirun $E/bin/mpiexec 2>/dev/null | head -1)
mkdir -p ~/dcv/gard_mpi2 && cd ~/dcv/gard_mpi2 && cp ../gard_input.fasta .
echo "launcher: $MPI   env: $E   ranks: ${NP:-17}"
$MPI -np ${NP:-17} $E/bin/HYPHYMPI gard --alignment gard_input.fasta --type nucleotide --rv Gamma --mode Faster --output gard_result.json --output-lf gard_fit.bf 2>&1 | tee gard.log
