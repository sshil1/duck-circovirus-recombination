#!/usr/bin/env python3
"""STEP 5a - build the analysis set for recombination / tree and align it.
Set A = full-length (>=1960 bp) genomes, EXCLUDING the lone DuCV-3 record (OP432310), which is ~73% identical to everything else:
including it would create ~1.7 kb of alignment gaps and dominate every distance-based test. It is analysed separately.
The 54 near-complete records (~1.88 kb) are held out of Set A for now (missing ~110 bp each); step 5b can add them back.
Outputs: setA_unaligned.fasta, tree_aln.fasta (MAFFT, columns with <95% occupancy removed), setA_ids.txt, setA_summary.txt
Usage: python3 step5_sets.py [threads]
"""
import csv, subprocess, sys
from collections import Counter
from Bio import SeqIO
T = sys.argv[1] if len(sys.argv) > 1 else "4"
DUCV3 = {"OP432310.1"}
meta = {r["accession"]: r for r in csv.DictReader(open("dcv_curated.tsv"), delimiter="\t")}
keep = []
for r in SeqIO.parse("dcv_oriented.fasta", "fasta"):
    acc = r.id.split("|")[0]; m = meta[acc]
    if m["completeness"] == "full" and acc not in DUCV3: keep.append(r)
SeqIO.write(keep, "setA_unaligned.fasta", "fasta")
open("setA_ids.txt", "w").write("\n".join(r.id for r in keep))
with open("setA_aln_raw.fasta", "w") as o: subprocess.run(["mafft", "--auto", "--thread", T, "--quiet", "setA_unaligned.fasta"], stdout=o, check=True)
aln = [(r.id, str(r.seq).upper()) for r in SeqIO.parse("setA_aln_raw.fasta", "fasta")]   # MAFFT lowercases: .upper() matters
L = len(aln[0][1]); n = len(aln)
occ = [sum(s[j] != "-" for _, s in aln) / n for j in range(L)]
cols = [j for j in range(L) if occ[j] >= 0.95]
with open("tree_aln.fasta", "w") as o:
    for i, s in aln: o.write(f">{i}\n{''.join(s[j] for j in cols)}\n")
var = sum(len({s[j] for _, s in aln} - {"-", "N"}) > 1 for j in cols)
txt = (f"Set A: {n} genomes   alignment {L} cols -> {len(cols)} after >=95% occupancy trim   variable sites {var}\n"
       f"countries: {dict(Counter(meta[i.split('|')[0]]['country'] or 'NA' for i,_ in aln).most_common())}\n")
open("setA_summary.txt", "w").write(txt); print(txt)
