#!/usr/bin/env bash
set -uo pipefail
H=/opt/conda/envs/hyphy/bin/hyphy
cd ~/dcv/sel
run() {
  P=$1
  for M in fubar fel busted meme; do
    echo "== $P $M $(date +%T)"
    EXTRA=""
    [ "$M" = busted ] && EXTRA="--srv Yes"
    [ "$M" = fel ] && EXTRA="--pvalue 0.1"
    [ "$M" = meme ] && EXTRA="--pvalue 0.1"
    $H CPU=9 "ENV=TOLERATE_NUMERICAL_ERRORS=1;" $M --alignment $P.fasta --tree $P.treefile $EXTRA --output hyphy/$P.$M.json > hyphy/$P.$M.log 2>&1 \
      && echo "done $P $M $(date +%T)" || echo "FAILED $P $M $(date +%T)"
  done
}
export -f run; export H
printf "%s\n" repA repB capA capB | xargs -P 4 -I{} bash -c 'run {}'
echo ALL_DONE $(date +%T)
