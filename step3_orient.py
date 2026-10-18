#!/usr/bin/env python3
"""STEP 3 - put every genome on the same strand and the same start coordinate (circular genomes are deposited with arbitrary starts/strands).
Method: pick the reference = most common-length record; for each genome and each strand, vote on the circular offset using shared 14-mers; keep the strand with more votes.
Outputs: dcv_oriented.fasta, dcv_orient_report.txt (flags genomes whose vote is weak = possibly divergent, e.g. DuCV-3, NOT errors).
"""
import re, sys
from collections import Counter, defaultdict
from Bio import SeqIO
K = 14
recs = [(r.id, r.description, str(r.seq).upper()) for r in SeqIO.parse("dcv_core.fasta", "fasta")]
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
L = Counter(len(s) for _, _, s in recs).most_common(1)[0][0]
ref_id, _, ref = next(x for x in recs if len(x[2]) == L)
idx = defaultdict(list)
for i in range(len(ref)): idx[(ref + ref[:K])[i:i+K]].append(i)
def vote(s):
    ss = s + s[:K - 1]; c = Counter()
    for i in range(len(s)):
        for j in idx.get(ss[i:i+K], ()): c[(i - j) % len(s)] += 1
    return c.most_common(1)[0] if c else (0, 0)
out = open("dcv_oriented.fasta", "w"); rep = []; flipped = rotated = weak = 0
for rid, desc, s in recs:
    off_f, v_f = vote(s)[0], vote(s)[1]; off_r, v_r = vote(rc(s))[0], vote(rc(s))[1]
    strand, off, v, seq = ("+", off_f, v_f, s) if v_f >= v_r else ("-", off_r, v_r, rc(s))
    seq = seq[off:] + seq[:off]
    flipped += strand == "-"; rotated += off != 0; weak += v < 100
    rep.append((rid, strand, off, v, len(s)))
    out.write(f">{desc}\n{seq}\n")
out.close()
with open("dcv_orient_report.txt", "w") as f:
    f.write(f"reference {ref_id} length {L}\nn={len(recs)}  reverse-complemented={flipped}  rotated={rotated}  weak-vote(<100 shared 14-mers)={weak}\n")
    f.write("weak-vote records (inspect; may be genuinely divergent lineage):\n")
    for r in sorted(rep, key=lambda x: x[3])[:30]: f.write("  %s strand=%s offset=%d votes=%d len=%d\n" % r)
print(open("dcv_orient_report.txt").read())
