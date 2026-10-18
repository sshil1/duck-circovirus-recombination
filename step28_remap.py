import csv
from Bio import SeqIO
raw={r.id:str(r.seq).upper() for r in SeqIO.parse("setA_aln_raw.fasta","fasta")}
tre={r.id:str(r.seq).upper() for r in SeqIO.parse("tree_aln.fasta","fasta")}
ids=list(raw); n=len(ids); L=len(raw[ids[0]]); common=[i for i in tre if i in raw][:30]; cols=None
for ch in ("-","-N"):
    c=[j for j in range(L) if sum(raw[i][j] not in ch for i in ids)/n>=0.95]
    if len(c)==len(next(iter(tre.values()))) and all("".join(raw[i][j] for j in c)==tre[i] for i in common): cols=c; break
assert cols, "could not reproduce tree_aln column set - paste this output"
ref=next(i for i in raw if i.startswith("HM162346")); coord=[]; k=0
for chh in raw[ref]:
    if chh!="-": k+=1
    coord.append(k)
mp=lambda c: coord[cols[c-1]]
print("kept",len(cols),"of",L,"columns | reference",ref.split("|")[0],"| column->HM162346.1:",{c:mp(c) for c in (1,578,680,1116,1640,1980)})
REP,CAP=(48,927),(1156,1930)
def region(p):
    if REP[0]<=p<=REP[1]: return "rep"
    if CAP[0]<=p<=CAP[1]: return "cap"
    if REP[1]<p<CAP[0]: return "intergenic (rep-cap 3' ends)"
    return "intergenic (origin region)"
rows=list(csv.reader(open("events_table.tsv"),delimiter="\t")); hdr=rows[0]
i1,i2,r1,r2=[hdr.index(x) for x in ("breakpoint1","breakpoint2","bp1_region","bp2_region")]
hdr=[("3SEQ_breakpoints_alncol" if h=="3SEQ_breakpoints" else h) for h in hdr]
out=[hdr+["breakpoint1_alncol","breakpoint2_alncol"]]
for r in rows[1:]:
    new=list(r); orig=[r[i1],r[i2]]
    for ci,ri in ((i1,r1),(i2,r2)):
        a,b=map(int,r[ci].split("-")); A,B=mp(a),mp(b); new[ci]=f"{A}-{B}"
        if r[ri]!="origin junction": new[ri]=region((A+B)/2)
    out.append(new+orig)
with open("events_table_HM.tsv","w") as f:
    for r in out: f.write("\t".join(r)+"\n")
for r in out: print("\t".join(r[x] for x in (0,3,10,11,12,13)))
src=open("step27_fig.py").read()
for old,new in (("GEN=1995","GEN=1993"),("((928,1158),(1933,GEN),(1,48))","((928,1155),(1931,GEN),(1,47))"),
    ("FancyArrow(49,0,878,0,","FancyArrow(48,0,879,0,"),("FancyArrow(1932,0,-773,0,","FancyArrow(1930,0,-774,0,"),
    ('open("events_table.tsv")','open("events_table_HM.tsv")'),("GQ423742.1 numbering","HM162346.1 numbering"),
    ("fig_breakpoint_map_v2","fig_breakpoint_map_v3"),("marked at 1-25","at the left edge"),
    ("(3SEQ places it at 1951-1980)","(3SEQ places it in the last 30 alignment columns)")):
    assert old in src, "patch target missing: "+old
    src=src.replace(old,new)
open("step27b_fig.py","w").write(src)
