import random
from Bio import SeqIO
S={r.id.split("|")[0]:r for r in SeqIO.parse("tree_aln.fasta","fasta")}
rc=[l.split("\t") for l in open("recomb_candidates.tsv")][1:]
flag=[r[0] for r in rc if float(r[3])>0.5]
parents="GQ423742.1 DQ166837.1 FJ554673.1 GQ423745.1 MN928809.1 MN928805.1 MN928807.1 MT646346.1 GQ423743.1 OR134499.1 OR134469.1 MN928808.1 PX471097.1 KC851822.1".split()
gt={l.split("\t")[0]:l.split("\t")[3] for l in list(open("genotype_by_segment.tsv"))[1:]}
random.seed(7)
bg=[]
for g in sorted(set(gt.values())):
    pool=[a for a in gt if gt[a]==g and a not in flag and a not in parents]
    bg+=random.sample(pool,min(8,len(pool)))
ids=[i for i in dict.fromkeys(flag+parents+bg) if i in S]
for i in ids: S[i].id=i; S[i].description=""
SeqIO.write([S[i] for i in ids],"rec3_subset.fasta","fasta")
print(len(flag),"flagged;",len(ids),"sequences written to rec3_subset.fasta")
