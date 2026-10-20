#!/usr/bin/env bash
# AU tests for sub-groups of group 1, per arc, reusing g1_out/arcI.* and arcO.* from g1_monophyly.sh.
# Sub-groups: Henan6 (clusters A+B), Fujian4, ClusterC (EU344803 + EU499311).
# For each arc x sub-group: best tree with the sub-group forced monophyletic vs the unconstrained ML tree, AU test.
# Run in ~/dcv:
#   tmux new -s sub "bash -c 'eval \"\$(conda shell.bash hook)\"; conda activate cavphylo; bash g1_subgroups.sh 2>&1 | tee g1_subgroups.log; exec bash'"
set -uo pipefail
OUT=g1_out; T=${T:-16}
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
declare -A SG
SG[Henan6]="MN928792.1 MN928802.1 MN928803.1 MN928804.1 MN928806.1 MN928811.1"
SG[Fujian4]="EU344802.1 EU344803.1 EU344805.1 EU499311.1"
SG[ClusterC]="EU344803.1 EU499311.1"
SUM=$OUT/subgroups_summary.tsv; printf "arc\tsubgroup\tn\tdeltaL_constrained\tp_AU_constrained\tp_AU_unconstrained\n" > $SUM
for ARC in arcI arcO; do
  BEST=$(grep -m1 "Best-fit model" $OUT/$ARC.iqtree | awk '{print $NF}')
  for G in Henan6 Fujian4 ClusterC; do
    TAG=${ARC}_$G; echo; echo "=== $ARC / $G (model $BEST) ==="
    python3 - $OUT/$ARC.treefile $OUT/$TAG.constraint.tre "${SG[$G]}" <<'P'
import sys
from Bio import Phylo
tre,outf,acc=sys.argv[1],sys.argv[2],sys.argv[3].split()
t=Phylo.read(tre,"newick"); names=[c.name for c in t.get_terminals()]
grp=[n for n in names if any(n.startswith(a) for a in acc)]
assert len(grp)==len(acc),(len(grp),len(acc))
rest=[n for n in names if n not in grp]
open(outf,"w").write("(("+",".join(grp)+"),"+",".join(rest)+");\n")
P
    [ $? -eq 0 ] || { echo "constraint failed for $TAG"; continue; }
    $IQ -s $OUT/$ARC.fasta -m "$BEST" -g $OUT/$TAG.constraint.tre -T $T --prefix $OUT/${TAG}_con -redo --keep-ident -quiet
    cat $OUT/$ARC.treefile $OUT/${TAG}_con.treefile > $OUT/${TAG}_both.trees
    $IQ -s $OUT/$ARC.fasta -m "$BEST" -z $OUT/${TAG}_both.trees -zb 10000 -au -T $T --prefix $OUT/${TAG}_au -redo --keep-ident -quiet
    sed -n '/^Tree *logL/,/^$/p' $OUT/${TAG}_au.iqtree
    ROWS=$(sed -n '/^Tree *logL/,/^$/p' $OUT/${TAG}_au.iqtree | awk '$1==1||$1==2')
    D2=$(echo "$ROWS" | awk '$1==2{print $3}'); P2=$(echo "$ROWS" | awk '$1==2{print $(NF-1)}'); P1=$(echo "$ROWS" | awk '$1==1{print $(NF-1)}')
    N=$(echo ${SG[$G]} | wc -w); printf "%s\t%s\t%s\t%s\t%s\t%s\n" $ARC $G $N "$D2" "$P2" "$P1" >> $SUM
  done
done
echo; echo "== SUMMARY (p_AU_constrained < 0.05 = forcing the subgroup together is rejected) =="; column -t -s$'\t' $SUM
echo Done
