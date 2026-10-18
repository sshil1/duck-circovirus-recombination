import csv
from Bio import SeqIO
S={r.id:str(r.seq).upper() for r in SeqIO.parse("rec3_subset.fasta","fasta")}
rec={r["C_ACCNUM"]:r for r in csv.DictReader(open("seq3/3s.rec.csv"))}
rc=[l.rstrip("\n").split("\t") for l in open("recomb_candidates.tsv")][1:]
flag=[r[0] for r in rc if float(r[3])>0.5]
L=len(next(iter(S.values())))
print("child\tparents(P,Q)\tn_inf\tbest_single_parent\tbest_mosaic\tgain\tinside_segment\tinside_parent\t3SEQ_breakpoints")
for c in flag:
    r=rec[c]; P,Q=r["P_ACCNUM"],r["Q_ACCNUM"]; lab=[]
    for i in range(L):
        p,q,x=S[P][i],S[Q][i],S[c][i]
        if p==q or p not in "ACGT" or q not in "ACGT": continue
        if x==p: lab.append((i+1,0))
        elif x==q: lab.append((i+1,1))
    n=len(lab)
    if n<10: print(c,P+","+Q,n,"too few informative sites",sep="\t"); continue
    p0=[0]; p1=[0]
    for _,l in lab: p0.append(p0[-1]+(l==0)); p1.append(p1[-1]+(l==1))
    t0,t1=p0[-1],p1[-1]; base=max(t0,t1); best=(base,None,None,None)
    for i in range(n+1):
        for j in range(i+1,n+1):
            in0=p0[j]-p0[i]; in1=p1[j]-p1[i]
            for sc,who in ((t0-in0+in1,Q),(t1-in1+in0,P)):
                if sc>best[0]: best=(sc,i,j,who)
    sc,i,j,who=best
    seg=f"{lab[i][0]}-{lab[j-1][0]}" if i is not None else "-"
    print(c,P+","+Q,n,base,sc,sc-base,seg,who or "-",r["breakpoints"],sep="\t")
