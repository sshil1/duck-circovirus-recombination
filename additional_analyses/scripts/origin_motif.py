#!/usr/bin/env python3
"""Where is the replication-origin nonanucleotide / stem-loop of DuCV (HM162346.1 numbering),
and how close are the recombination breakpoints to it?

Reads the alignment (default tree_aln.fasta), takes HM162346.1, strips gaps, then
  1. lists every NANTATTAC match (both strands) - the circovirus-family origin nonamer is a
     degenerate version of this; we do not assume which hit is functional,
  2. finds the strongest hairpins (stem >= 8 bp, loop 3-12 nt, >=1 GC pair weighting) and whether
     one carries a nonamer in its loop,
  3. gives the circular distance of each breakpoint interval to the nearest nonamer / hairpin,
  4. compares with random breakpoint positions (20,000 draws, same number of breakpoints).
Usage: python3 origin_motif.py [tree_aln.fasta]
"""
import re, sys, numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
fn = sys.argv[1] if len(sys.argv) > 1 else "tree_aln.fasta"
ref = None
for r in SeqIO.parse(fn, "fasta"):
    if r.id.startswith("HM162346.1"):
        ref = str(r.seq).upper().replace("-", ""); break
assert ref, "HM162346.1 not found"
L = len(ref); print(f"HM162346.1: {L} nt (circular)")
ext = ref + ref[:30]          # allow motifs across the origin
rc = str(Seq(ref).reverse_complement()); extrc = rc + rc[:30]

hits = []
for m in re.finditer(r"(?=(.A.TATTAC))", ext):
    p = m.start() + 1
    if p <= L: hits.append((p, "+", ext[m.start():m.start()+9]))
for m in re.finditer(r"(?=(.A.TATTAC))", extrc):
    p = m.start() + 1
    if p <= L: hits.append((L - p - 7, "-", extrc[m.start():m.start()+9]))   # leftmost position on + strand
print("\nNANTATTAC matches (position = leftmost nt on the + strand, 1-based):")
for p, s, q in sorted(hits): print(f"  {p:5d}  strand {s}  {q}")

# hairpins: stems >= 6 bp (G-T wobble allowed) with loop 3-12, scored 3 per G/C and 2 per A/T in the stem
pairs = {("A","T"),("T","A"),("G","C"),("C","G"),("G","T"),("T","G")}
def hairpins(s, minstem=6):
    n = len(s); out = []
    for i in range(n):
        for loop in range(3, 13):
            k = 0
            while k < 25:
                a, b = i - k, i + 1 + loop + k      # a = left stem nt, b = right stem nt (5'->3' pairs a with b)
                if a < 0 or b >= n or (s[a], s[b]) not in pairs: break
                k += 1
            if k >= minstem:
                stem = s[i - k + 1:i + 1]
                sc = sum(3 if c in "GC" else 2 for c in stem)   # wobble G-U/T pairs counted as 2
                out.append((sc, i - k + 1 + 1, i + 1 + loop + k, k, loop, s[i + 1:i + 1 + loop]))
    return sorted(out, reverse=True)
hp = hairpins(ref)
print("\nstrongest perfect hairpins (score, start..end, stem bp, loop nt, loop seq):")
seen = []; hp_centres = []
for sc, a, b, k, loop, ls in hp:
    if any(abs(a - x) < 10 for x in seen): continue
    seen.append(a); hp_centres.append((a + b) // 2)
    print(f"  score {sc:3d}  {a:5d}..{b:5d}  stem {k:2d} bp  loop {loop:2d}  {ls}")
    if len(seen) >= 6: break

# breakpoint intervals (HM162346.1 numbering) from events_table_HM.tsv, one row per event cluster
events = {
 "A MN928802/803/804/806": ((1, 25), (579, 675)),
 "B MN928792/MN928811":    ((1, 25), (579, 726)),
 "C EU344803/EU499311":    ((60, 72), (579, 682)),
 "EU344802":               ((1, 25), (579, 682)),
 "EU344805":               ((46, 54), (672, 682)),
 "OR134497":               ((1110, 1125), (1649, 1723)),
 "OR134498":               ((871, 975), (1655, 1664)),
 "OK094644":               ((198, 222), (582, 672)),
 "KP943594":               ((180, 192), (1125, 1138)),
}
motif_pos = np.array(sorted({p for p, s, q in hits}))
top_hp = np.array(sorted(hp_centres))
def cdist(x, pts):
    if len(pts) == 0: return float("nan")
    d = np.abs(pts - x); return np.minimum(d, L - d).min()
mids = []
print("\nbreakpoint interval midpoints and circular distance (nt) to nearest nonamer / strong hairpin centre:")
for name, (b1, b2) in events.items():
    for b in (b1, b2):
        m = (b[0] + b[1]) / 2; mids.append(m)
        print(f"  {name:24s} {b[0]:5d}-{b[1]:<5d} mid {m:7.1f}  nonamer {cdist(m, motif_pos):5.0f}   hairpin {cdist(m, top_hp):5.0f}")
mids = np.array(mids)
rng = np.random.default_rng(7)
obs = np.mean([cdist(m, motif_pos) for m in mids])
null = np.array([np.mean([cdist(x, motif_pos) for x in rng.uniform(0, L, len(mids))]) for _ in range(20000)])
print(f"\nmean distance of the {len(mids)} breakpoints to the nearest nonamer match: {obs:.1f} nt; "
      f"random placement: {null.mean():.1f} nt; P(random <= observed) = {(np.sum(null <= obs)+1)/20001:.4f}")
print("NB: with several NANTATTAC matches in a 2-kb genome this is a weak test; use it to locate the origin, not to infer mechanism.")
