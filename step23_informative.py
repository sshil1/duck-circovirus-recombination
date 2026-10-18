from Bio import SeqIO
from collections import Counter
S={r.id:str(r.seq).upper() for r in SeqIO.parse("rec3_subset.fasta","fasta")}
R6="MN928811.1 MN928792.1 MN928806.1 MN928803.1 MN928804.1 MN928802.1".split()
A=["GQ423745.1","GQ423743.1"]; D=["MN928807.1","MN928808.1","MN928809.1"]
L=len(next(iter(S.values()))); rows=[]
for i in range(L):
    a={S[x][i] for x in A}; d={S[x][i] for x in D}
    if len(a)!=1 or len(d)!=1: continue
    ab,db=a.pop(),d.pop()
    if ab==db or ab not in "ACGT" or db not in "ACGT": continue
    rb=[S[x][i] for x in R6 if S[x][i] in "ACGT"]
    if len(rb)<5: continue
    top=max(set(rb),key=rb.count)
    if rb.count(top)<5: continue
    rows.append((i+1,"A" if top==ab else "D" if top==db else "x"))
print("informative sites (A group != D group):",len(rows))
for n,a,b in [("1-577",1,577),("578-724",578,724),("725-1950",725,1950),("1951-1980",1951,1980)]:
    print(f"{n:10s}",dict(Counter(l for p,l in rows if a<=p<=b)))
print("\nsites 450-900:"," ".join(f"{p}{l}" for p,l in rows if 450<=p<=900))
