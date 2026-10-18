#!/usr/bin/env bash
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh
E=/opt/conda/envs/hyphy
if [ ! -x $E/bin/mpirun ] && [ ! -x $E/bin/mpiexec ]; then
  echo "no MPI launcher in hyphy env -> building matched env hyphympi"
  conda create -y -n hyphympi -c conda-forge -c bioconda hyphy openmpi
  E=/opt/conda/envs/hyphympi
fi
MPI=$(ls $E/bin/mpirun $E/bin/mpiexec 2>/dev/null | head -1)
mkdir -p ~/dcv/gard_mpi && cd ~/dcv/gard_mpi && cp ../gard_input.fasta .
N=${NP:-17}
echo "launcher: $MPI   ranks: $N (1 master + $((N-1)) workers)"
$MPI -np $N $E/bin/HYPHYMPI gard --alignment gard_input.fasta --type nucleotide --rv GDD --rate-classes 3 --mode Faster --output gard_result.json --output-lf gard_fit.bf 2>&1 | tee gard.log
