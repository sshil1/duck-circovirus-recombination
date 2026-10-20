#!/usr/bin/env bash
# Is "group 1" (10 genomes sharing the right breakpoint) one lineage or several?
# Builds ML trees of the two non-recombinant arcs, tests monophyly of group 1 and of its
# Fujian-4 / Henan-6 sub-groups, and runs an AU test of the constrained vs unconstrained tree.
# Run:  tmux new -s g1 'bash g1_monophyly.sh tree_aln.fasta 2>&1 | tee g1.log'
# Needs: python3+biopython, iqtree3 (or iqtree2) on PATH (server: conda activate cavphylo).
set -euo pipefail
ALN=${1:-tree_aln.fasta}; T=${T:-16}; OUT=g1_out; mkdir -p $OUT
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
# arc I (between the two breakpoints of group 1, alignment cols) and arc O (outside both)
ARC_I=${ARC_I:-80-570}; ARC_O=${ARC_O:-730-1940}
G1="EU344802.1 EU344803.1 EU344805.1 EU499311.1 MN928792.1 MN928802.1 MN928803.1 MN928804.1 MN928806.1 MN928811.1"
FJ="EU344802.1 EU344803.1 EU344805.1 EU499311.1"; HN="MN928792.1 MN928802.1 MN928803.1 MN928804.1 MN928806.1 MN928811.1"
python3 - "$ALN" "$ARC_I" "$ARC_O" "$OUT" <<'P'
import sys
from Bio import SeqIO
aln,ai,ao,out=sys.argv[1:5]
recs=list(SeqIO.parse(aln,"fasta"))
for tag,rg in (("arcI",ai),("arcO",ao)):
    a,b=map(int,rg.split("-"))
    with open(f"{out}/{tag}.fasta","w") as f:
        for r in recs: f.write(f">{r.id}\n{str(r.seq)[a-1:b]}\n")
    print(tag,rg,b-a+1,"cols")
P
for ARC in arcI arcO; do
  $IQ -s $OUT/$ARC.fasta -m ${MODEL:-MFP} -B ${BOOT:-1000} -T $T --prefix $OUT/$ARC -redo --keep-ident -quiet
done
python3 - "$OUT" "$G1" "$FJ" "$HN" <<'P'
import sys, re
from Bio import Phylo
out=sys.argv[1]; sets={"group1(10)":sys.argv[2].split(),"Fujian(4)":sys.argv[3].split(),"Henan(6)":sys.argv[4].split()}
for arc in ("arcI","arcO"):
    t=Phylo.read(f"{out}/{arc}.treefile","newick"); t.root_at_midpoint()
    tips={c.name:c for c in t.get_terminals()}
    print(f"\n== {arc}: {len(tips)} tips ==")
    for lab,acc in sets.items():
        names=[n for n in tips if n.split("|")[0].split("_")[0] in acc or any(n.startswith(a) for a in acc)]
        if len(names)!=len(acc): print(lab,"WARNING matched",len(names),"of",len(acc),names)
        m=t.common_ancestor([tips[n] for n in names]); n_tips=len(m.get_terminals())
        sup=m.confidence
        print(f"{lab:11s} MRCA clade has {n_tips} tips ({'MONOPHYLETIC' if n_tips==len(names) else 'NOT monophyletic: %d extra tips'%(n_tips-len(names))}); support={sup}")
P
# constrained trees + AU test, per arc
for ARC in arcI arcO; do
python3 - $OUT $ARC "$G1" <<'P'
import sys
from Bio import Phylo
out,arc,g1=sys.argv[1],sys.argv[2],sys.argv[3].split()
t=Phylo.read(f"{out}/{arc}.treefile","newick")
names=[c.name for c in t.get_terminals()]
grp=[n for n in names if any(n.startswith(a) for a in g1)]
rest=[n for n in names if n not in grp]
open(f"{out}/{arc}.constraint.tre","w").write("(("+",".join(grp)+"),"+",".join(rest)+");\n")
P
  BEST=$(grep -m1 "Best-fit model" $OUT/$ARC.iqtree | sed 's/.*: *\([^ ]*\) chosen.*/\1/'); [ -n "$BEST" ] || BEST=${MODEL:-GTR+G}
  $IQ -s $OUT/$ARC.fasta -m "$BEST" -g $OUT/$ARC.constraint.tre -T $T --prefix $OUT/${ARC}_con -redo --keep-ident -quiet
  cat $OUT/$ARC.treefile $OUT/${ARC}_con.treefile > $OUT/${ARC}_both.trees
  $IQ -s $OUT/$ARC.fasta -m "$BEST" -z $OUT/${ARC}_both.trees -zb 10000 -au -T $T --prefix $OUT/${ARC}_au -redo --keep-ident -quiet
  echo; echo "== AU test, $ARC (tree 1 = unconstrained ML, tree 2 = group 1 forced monophyletic) =="
  sed -n '/^USER TREES/,/^TIME STAMP/p' $OUT/${ARC}_au.iqtree | head -20
done
echo; echo "Done. Key output: g1_out/arcI.treefile, g1_out/arcO.treefile and the AU tables above."
