#!/usr/bin/env python3
"""STEP 4 - align (MAFFT), compute pairwise identity, and test the taxonomy question: how many distinct groups exist at the ICTV 80% species threshold, and what is each?
Outputs: dcv_aln.fasta, dcv_identity.tsv, dcv_clusters.txt, dcv_nn.tsv (nearest-neighbour identity per genome)
Usage: python3 step4_cluster.py [threads]
"""
import subprocess, sys, numpy as np
from collections import Counter
from Bio import SeqIO
T = sys.argv[1] if len(sys.argv) > 1 else "4"
with open("dcv_aln.fasta", "w") as o: subprocess.run(["mafft", "--auto", "--thread", T, "--quiet", "dcv_oriented.fasta"], stdout=o, check=True)
recs = [(r.id, str(r.seq).upper()) for r in SeqIO.parse("dcv_aln.fasta", "fasta")]
ids = [r[0] for r in recs]; A = np.array([list(s) for _, s in recs]); n = len(ids)
code = {c: i for i, c in enumerate("ACGT")}; M = np.full(A.shape, -1, dtype=np.int8)
for c, i in code.items(): M[A == c] = i
ident = np.zeros((n, n))
for i in range(n):
    valid = (M[i] >= 0) & (M >= 0); same = (M == M[i]) & valid
    ident[i] = same.sum(1) / np.maximum(valid.sum(1), 1)
np.fill_diagonal(ident, 1.0)
with open("dcv_identity.tsv", "w") as f:
    f.write("id\t" + "\t".join(ids) + "\n"); [f.write(ids[i] + "\t" + "\t".join(f"{x:.4f}" for x in ident[i]) + "\n") for i in range(n)]
def comps(th):
    par = list(range(n))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(n):
        for j in np.where(ident[i] >= th)[0]:
            par[find(i)] = find(j)
    g = {}; [g.setdefault(find(i), []).append(i) for i in range(n)]
    return sorted(g.values(), key=len, reverse=True)
with open("dcv_nn.tsv", "w") as f:
    f.write("id\tnearest_identity\tnearest_id\n")
    for i in range(n):
        r = ident[i].copy(); r[i] = -1; j = int(r.argmax()); f.write(f"{ids[i]}\t{r[j]:.4f}\t{ids[j]}\n")
out = []; P = out.append
P(f"n={n}  alignment columns={A.shape[1]}  median pairwise identity={np.median(ident[np.triu_indices(n,1)]):.4f}  min={ident[np.triu_indices(n,1)].min():.4f}")
for th in (0.80, 0.85, 0.90, 0.95):
    cs = comps(th); P(f"\nsingle-linkage components at {int(th*100)}% identity: {len(cs)}  sizes: {[len(c) for c in cs][:15]}")
    if th in (0.80, 0.85):
        for k, c in enumerate(cs[:6]):
            P(f"  comp{k+1} n={len(c)}  e.g. {', '.join(ids[i].split('|')[0].split()[0] for i in c[:5])}  hosts={dict(Counter(ids[i].split('|')[2] for i in c if ids[i].count('|')>=2))}")
cs = comps(0.80); main = set(cs[0])
others = [i for i in range(n) if i not in main]
P(f"\ngenomes outside the largest 80%-component: {len(others)}")
for i in others[:40]:
    best = max(ident[i, j] for j in main); P(f"  {ids[i]}  best identity to main component {best:.3f}")
open("dcv_clusters.txt", "w").write("\n".join(out)); print("\n".join(out))
