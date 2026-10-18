import sys
from Bio import SeqIO, Entrez
ref="GQ423742.1"; cols=[578,680,1116,1640]
aln=list(SeqIO.parse("rec3_subset.fasta","fasta"))
s=str([x for x in aln if x.id==ref][0].seq).upper()
m={}; p=0
for i,ch in enumerate(s,1):
    if ch!="-": p+=1
    m[i]=p if ch!="-" else None
ung=s.replace("-","")
print(f"{ref}: ungapped length {len(ung)}; first non-gap column {next(i for i,c in enumerate(s,1) if c!='-')}")
for c in cols:
    lo,hi=max(1,c-10),c+10
    g=sum(str(x.seq)[lo-1:hi].count("-") for x in aln)/(len(aln)*(hi-lo+1))
    print(f"column {c}: position in {ref} = {m[c]} | codon position {(c-1)%3+1} | gap fraction +/-10 cols = {g:.3f}")
Entrez.email=sys.argv[1]
try:
    rec=SeqIO.read(Entrez.efetch(db="nucleotide",id=ref,rettype="gb",retmode="text"),"genbank")
    gs=str(rec.seq).upper(); off=gs.find(ung[:40])
    print("\nGenBank:",rec.description,"| record length",len(gs))
    print("alignment region starts at GenBank position", off+1 if off>=0 else "NOT FOUND (trimmed/reverse-complemented?)")
    for f in rec.features:
        if f.type in ("gene","CDS","mat_peptide","misc_feature","5'UTR","3'UTR"):
            q=f.qualifiers; print(f.type, f.location, q.get("gene",[""])[0], q.get("product",[""])[0])
    if off>=0:
        print("\nbreakpoint columns as GenBank positions:",{c:(m[c]+off if m[c] else None) for c in cols})
except Exception as e:
    print("\nNCBI fetch failed:",e)
