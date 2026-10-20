#!/usr/bin/env bash
# Selection sensitivity: rerun FUBAR/FEL/MEME/BUSTED on repA/repB/capA/capB after removing ALL 59 flagged
# genomes (recomb_candidates.tsv flag==FLAG; the original run removed only the 14 with score_main>0.5),
# then list sites supported by >=2 methods and compare with the seven candidate sites from the first run.
# Run from ~/dcv:   tmux new -s sel 'bash sel_noflag.sh 2>&1 | tee sel_noflag.log'
set -euo pipefail
cd ~/dcv; D=sel_noflag; mkdir -p $D/hyphy; 
source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
python3 - <<'P'
from Bio import SeqIO
rows=[l.rstrip("\n").split("\t") for l in open("recomb_candidates.tsv")]
h=rows[0]; fi=h.index("flag"); bad={r[0].split(".")[0] for r in rows[1:] if r[fi]=="FLAG"}
print("flagged genomes removed:",len(bad))
REF="HM162346.1"
def load(fn): return {r.id.split("|")[0]:str(r.seq).upper() for r in SeqIO.parse(fn,"fasta")}
def split(fn,nt,names):
    S=load(fn); ref=S[REF]; c=0
    for j,ch in enumerate(ref):
        if ch!="-": c+=1
        if c==nt: cut=j+1; break
    cut=((cut+2)//3)*3
    for nm,(a,b) in zip(names,[(0,cut),(cut,len(ref))]):
        seqs={k:v[a:b] for k,v in S.items() if k.split(".")[0] not in bad}
        seen={}
        for k in sorted(seqs,key=lambda x:(x!=REF,x)): seen.setdefault(seqs[k],k)
        keep=sorted(seen.values())
        with open(f"sel_noflag/{nm}.fasta","w") as o:
            for k in keep: o.write(f">{k.replace('.','_')}\n{seqs[k]}\n")
        print(nm,(b-a)//3,"codons;",len(seqs),"->",len(keep),"unique seqs")
split("rep_codon_aln.fasta",630,["repA","repB"]); split("cap_codon_aln.fasta",234,["capA","capB"])
P
IQ=$(command -v iqtree3 || command -v iqtree2 || (conda activate cavphylo; command -v iqtree3))
for P in repA repB capA capB; do $IQ -s $D/$P.fasta -m GTR+F+R4 -T 8 --prefix $D/$P -redo -quiet; cp $D/$P.treefile $D/$P.tre; done
H=/opt/conda/envs/hyphy/bin/hyphy
run(){ P=$1; M=$2; X=""; [ $M = busted ] && X="--srv Yes"; { [ $M = fel ] || [ $M = meme ]; } && X="--pvalue 0.1"
  $H CPU=2 "ENV=TOLERATE_NUMERICAL_ERRORS=1;" $M --alignment $D/$P.fasta --tree $D/$P.treefile $X --output $D/hyphy/$P.$M.json > $D/hyphy/$P.$M.log 2>&1; }
export -f run; export H D
printf '%s\n' repA repB capA capB | while read P; do for M in fubar fel busted meme; do echo "$P $M"; done; done | xargs -P 4 -L1 bash -c 'run $0 $1'
sed 's#~/dcv/sel"#~/dcv/sel_noflag"#' step17_summary.py > $D/summary.py   # needs step17_summary.py in ~/dcv
python3 $D/summary.py | tee $D/summary.txt
echo; echo "Compare the listed 'sites supported by >=2 methods' with the seven candidate sites from the original run (step17 output)."
