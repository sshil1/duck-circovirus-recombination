import csv
from Bio import SeqIO
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, Patch
S={r.id:str(r.seq).upper() for r in SeqIO.parse("rec3_subset.fasta","fasta")}
rec={r["C_ACCNUM"]:r for r in csv.DictReader(open("seq3/3s.rec.csv"))}
rcr=[l.rstrip("\n").split("\t") for l in open("recomb_candidates.tsv")][1:]
meta={r[0]:(r[1],r[2]) for r in rcr}; flag=[r[0] for r in rcr if float(r[3])>0.5]
L=len(next(iter(S.values()))); GEN=1995
# author's grouping by breakpoint signature (edit if you disagree)
GROUP={}
for a in "MN928811.1 MN928792.1 MN928806.1 MN928803.1 MN928804.1 MN928802.1 EU344802.1 EU344805.1 EU499311.1 EU344803.1".split(): GROUP[a]="1 rep(~578-680)/origin junction"
for a in "OR134497.1 OR134498.1".split(): GROUP[a]="2 intergenic/cap"
GROUP["OK094644.1"]="3 within rep"; GROUP["KP943594.1"]="4 rep/intergenic (poor fit)"
def region(p):
    p=int(round(p))
    if 49<=p<=927: return "rep"
    if 1159<=p<=1932: return "cap"
    if 928<=p<=1158: return "intergenic (rep-cap 3' ends)"
    return "intergenic (origin region)"
rows=[]
for c in flag:
    r=rec[c]; P,Q=r["P_ACCNUM"],r["Q_ACCNUM"]; lab=[]
    for i in range(L):
        p,q,x=S[P][i],S[Q][i],S[c][i]
        if p==q or p not in "ACGT" or q not in "ACGT": continue
        if x==p: lab.append((i+1,0))
        elif x==q: lab.append((i+1,1))
    n=len(lab); p0=[0]; p1=[0]
    for _,l in lab: p0.append(p0[-1]+(l==0)); p1.append(p1[-1]+(l==1))
    t0,t1=p0[-1],p1[-1]; base=max(t0,t1); best=(base,None,None,None,None)
    for i in range(n+1):
        for j in range(i+1,n+1):
            in0=p0[j]-p0[i]; in1=p1[j]-p1[i]
            for sc,ins,outs in ((t0-in0+in1,Q,P),(t1-in1+in0,P,Q)):
                if sc>best[0]: best=(sc,i,j,ins,outs)
    sc,i,j,ins,outs=best
    bp1=(lab[i-1][0] if i>0 else 1, lab[i][0]); bp2=(lab[j-1][0], lab[j][0] if j<n else L)
    m1=sum(bp1)/2; m2=sum(bp2)/2
    rows.append(dict(acc=c,country=meta[c][0],year=meta[c][1],group=GROUP.get(c,"?"),inside=ins,outside=outs,n=n,
        fit=100*sc/n,single=100*base/n,gain=sc-base,bp1=bp1,bp2=bp2,m1=m1,m2=m2,tseq=r["breakpoints"]))
rows.sort(key=lambda d:(d["group"],d["acc"]))
with open("events_table.tsv","w") as f:
    f.write("accession\tcountry\tyear\tgroup\tinside_parent\toutside_parent\tn_informative\tmosaic_fit_%\tbest_single_%\tgain_sites\tbreakpoint1\tbp1_region\tbreakpoint2\tbp2_region\t3SEQ_breakpoints\n")
    for d in rows:
        f.write("\t".join(map(str,[d["acc"],d["country"],d["year"],d["group"],d["inside"],d["outside"],d["n"],f'{d["fit"]:.1f}',f'{d["single"]:.1f}',d["gain"],
            f'{d["bp1"][0]}-{d["bp1"][1]}',region(d["m1"]) if d["bp1"][0]>1 else "origin junction",
            f'{d["bp2"][0]}-{d["bp2"][1]}',region(d["m2"]) if d["bp2"][1]<L else "origin junction",d["tseq"]]))+"\n")
fig,ax=plt.subplots(figsize=(11,0.34*len(rows)+2.2))
OR,BL,GR,PK="#E69F00","#0072B2","#009E73","#CC79A7"
for a,b in ((928,1158),(1933,GEN),(1,48)): ax.add_patch(Rectangle((a,-0.6),b-a,len(rows)+1.4,color="#eeeeee",zorder=0))
ax.add_patch(FancyArrow(49,0,878,0,width=0.5,head_width=0.7,head_length=40,length_includes_head=True,color=GR))
ax.add_patch(FancyArrow(1932,0,-773,0,width=0.5,head_width=0.7,head_length=40,length_includes_head=True,color=PK))
ax.text(480,0,"rep",ha="center",va="center",color="white",fontsize=9); ax.text(1545,0,"cap",ha="center",va="center",color="white",fontsize=9)
yl=[]
for k,d in enumerate(rows,1):
    ax.add_patch(Rectangle((1,k-0.3),GEN-1,0.6,color=BL))
    ax.add_patch(Rectangle((d["m1"],k-0.3),d["m2"]-d["m1"],0.6,color=OR))
    yl.append(f'{d["acc"]}  (in: {d["inside"]} | out: {d["outside"]})')
ax.set_yticks(range(1,len(rows)+1)); ax.set_yticklabels(yl,fontsize=7); ax.set_ylim(len(rows)+0.8,-0.9); ax.set_xlim(1,GEN)
ax.set_xlabel("Genome position (nt, GQ423742.1 numbering; circular genome, grey = intergenic regions)")
ax.legend(handles=[Patch(color=OR,label="arc from 'inside' parent"),Patch(color=BL,label="remainder (wraps through origin) from other parent")],loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=8,frameon=False)
for s in ("top","right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig("fig_breakpoint_map.png",dpi=300); plt.savefig("fig_breakpoint_map.pdf")
print(open("events_table.tsv").read())
