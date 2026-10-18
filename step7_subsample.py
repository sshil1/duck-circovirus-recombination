#!/usr/bin/env python3
"""STEP 7a - block-aware, diversity-maximising subsample of Set A for GARD (breakpoint detection).
Why: GARD cost explodes with taxa count, and 444 genomes are not 444 independent observations (98 submission blocks).
Rule 1: at most CAP genomes per submission block (the most mutually divergent ones).
Rule 2: if still > N, farthest-point (max-min distance) sampling down to N, seeded from the genome nearest the median.
Caveat (state it in the paper): farthest-point sampling over-represents divergent lineages; it is chosen to MAXIMISE breakpoint-detection power, not to estimate frequencies.
Outputs: gard_input.fasta (sanitised ids = accession), gard_id_map.tsv, subsample_summary.txt
Usage: python3 step7_subsample.py [N=100] [CAP=2]
"""
import csv, sys, numpy as np
from collections import Counter
from Bio import SeqIO
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 2
meta = {r["accession"]: r for r in csv.DictReader(open("dcv_curated.tsv"), delimiter="\t")}
rows = [l.rstrip("\n").split("\t") for l in open("dcv_identity.tsv")]
ids_all = rows[0][1:]; D = 1 - np.array([[float(x) for x in r[1:]] for r in rows[1:]])
pos = {i: k for k, i in enumerate(ids_all)}
aln = {r.id: str(r.seq).upper() for r in SeqIO.parse("tree_aln.fasta", "fasta")}
setA = list(aln); acc = lambda i: i.split("|")[0]
blocks = {}
for i in setA: blocks.setdefault(meta[acc(i)]["block"] or acc(i), []).append(i)
chosen = []
for b, members in blocks.items():
    if len(members) <= CAP: chosen += members; continue
    sel = [members[0]]                                  # farthest-point inside the block
    while len(sel) < CAP:
        sel.append(max((m for m in members if m not in sel), key=lambda m: min(D[pos[m], pos[s]] for s in sel)))
    chosen += sel
n_after_cap = len(chosen)
if len(chosen) > N:
    idx = [pos[i] for i in chosen]; sub = D[np.ix_(idx, idx)]
    start = int(np.argmin(np.abs(sub.mean(1) - np.median(sub.mean(1)))))
    sel = [start]; mind = sub[start].copy()
    while len(sel) < N:
        k = int(np.argmax(mind)); sel.append(k); mind = np.minimum(mind, sub[k])
    chosen = [chosen[k] for k in sel]
chosen = sorted(chosen, key=acc)
with open("gard_input.fasta", "w") as o, open("gard_id_map.tsv", "w") as m:
    for i in chosen: o.write(f">{acc(i)}\n{aln[i]}\n"); m.write(f"{acc(i)}\t{i}\n")
d = D[np.ix_([pos[i] for i in chosen], [pos[i] for i in chosen])]
t = (f"Set A {len(setA)} genomes in {len(blocks)} submission blocks -> cap {CAP}/block gives {n_after_cap} -> final {len(chosen)}\n"
     f"countries: {dict(Counter(meta[acc(i)]['country'] or 'NA' for i in chosen).most_common())}\n"
     f"years: {dict(sorted(Counter(meta[acc(i)]['collection_year'] or 'NA' for i in chosen).items()))}\n"
     f"host: {dict(Counter(meta[acc(i)]['hostclass'] for i in chosen))}\n"
     f"blocks represented: {len({meta[acc(i)]['block'] for i in chosen})}   median pairwise distance in subsample {np.median(d[np.triu_indices(len(chosen),1)]):.4f} (Set A: {np.median(D[np.ix_([pos[i] for i in setA],[pos[i] for i in setA])][np.triu_indices(len(setA),1)]):.4f})\n")
open("subsample_summary.txt", "w").write(t); print(t)
