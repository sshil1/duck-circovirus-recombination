#!/usr/bin/env bash
# Re-uses g1_out/arcI.treefile and arcO.treefile (already built). Fixes two bugs of g1_monophyly.sh:
#  (1) monophyly is tested on the UNROOTED tree (root on a leaf outside the group, then ask whether the MRCA
#      clade contains exactly the group), not on a midpoint-rooted tree;
#  (2) the model name is read correctly for the constrained search and the AU test.
# Run in ~/dcv:  tmux new -s g1post "bash -c 'eval \"\$(conda shell.bash hook)\"; conda activate cavphylo; bash g1_post.sh 2>&1 | tee g1_post.log; exec bash'"
set -uo pipefail
OUT=g1_out; T=${T:-16}
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
G1="EU344802.1 EU344803.1 EU344805.1 EU499311.1 MN928792.1 MN928802.1 MN928803.1 MN928804.1 MN928806.1 MN928811.1"
FJ="EU344802.1 EU344803.1 EU344805.1 EU499311.1"; HN="MN928792.1 MN928802.1 MN928803.1 MN928804.1 MN928806.1 MN928811.1"
cat > $OUT/mono.py <<'P'
import sys
from Bio import Phylo
out=sys.argv[1]; sets={"group1(10)":sys.argv[2].split(),"Fujian(4)":sys.argv[3].split(),"Henan(6)":sys.argv[4].split()}
for arc in ("arcI","arcO"):
    t=Phylo.read(f"{out}/{arc}.treefile","newick")
    allt=[c.name for c in t.get_terminals()]
    print(f"\n== {arc}: {len(allt)} tips (unrooted test) ==")
    for lab,acc in sets.items():
        grp=[n for n in allt if any(n.startswith(a) for a in acc)]
        if len(grp)!=len(acc): print(lab,"WARNING matched",len(grp),"of",len(acc))
        outg=next(n for n in allt if n not in grp)
        t2=Phylo.read(f"{out}/{arc}.treefile","newick"); t2.root_with_outgroup(outg)
        tips={c.name:c for c in t2.get_terminals()}
        m=t2.common_ancestor([tips[n] for n in grp]); k=len(m.get_terminals())
        if k==len(grp): print(f"{lab:11s} MONOPHYLETIC (unrooted split); bootstrap/UFBoot support = {m.confidence}")
        else:
            extra=[c.name.split('|')[0] for c in m.get_terminals() if c.name not in grp][:8]
            print(f"{lab:11s} NOT monophyletic: smallest clade containing them has {k} tips ({k-len(grp)} outsiders, e.g. {', '.join(extra)})")
P
python3 $OUT/mono.py $OUT "$G1" "$FJ" "$HN"
for ARC in arcI arcO; do
  BEST=$(grep -m1 "Best-fit model" $OUT/$ARC.iqtree | awk '{print $NF}')
  echo; echo "== $ARC: model $BEST =="
  python3 - $OUT $ARC "$G1" <<'P'
import sys
from Bio import Phylo
out,arc,g1=sys.argv[1],sys.argv[2],sys.argv[3].split()
t=Phylo.read(f"{out}/{arc}.treefile","newick"); names=[c.name for c in t.get_terminals()]
grp=[n for n in names if any(n.startswith(a) for a in g1)]; rest=[n for n in names if n not in grp]
open(f"{out}/{arc}.constraint.tre","w").write("(("+",".join(grp)+"),"+",".join(rest)+");\n")
P
  $IQ -s $OUT/$ARC.fasta -m "$BEST" -g $OUT/$ARC.constraint.tre -T $T --prefix $OUT/${ARC}_con -redo --keep-ident -quiet
  cat $OUT/$ARC.treefile $OUT/${ARC}_con.treefile > $OUT/${ARC}_both.trees
  $IQ -s $OUT/$ARC.fasta -m "$BEST" -z $OUT/${ARC}_both.trees -zb 10000 -au -T $T --prefix $OUT/${ARC}_au -redo --keep-ident -quiet
  echo "== AU test, $ARC (tree 1 = unconstrained ML, tree 2 = group 1 forced monophyletic; p-AU < 0.05 rejects tree 2) =="
  sed -n '/^USER TREES/,/^TIME STAMP/p' $OUT/${ARC}_au.iqtree | head -20
done
echo; echo Done
