#!/usr/bin/env python3
"""Do breakpoints of independent DuCV recombination events cluster on the genome?
Unit = one event cluster (9 clusters, 2 breakpoints each, interval midpoints in alignment columns).
Statistic = largest number of the 18 breakpoints inside any window of width WIN on the circular genome.
Nulls (L = 1980 alignment columns, circular):
  uniform   each breakpoint uniform on the circle
  informative  each breakpoint drawn from the columns of the 370 biallelic informative sites
  rotation  each event's pair of breakpoints rotated as a rigid arc to a uniform position
Also reported at the level of the 4 breakpoint signatures (8 breakpoints).
Usage: python3 hotspot_test.py tree_aln.fasta
"""
import sys, numpy as np
from Bio import SeqIO
L = 1980; NSIM = 200000; rng = np.random.default_rng(42)
mid = lambda a, b: (a + b) / 2
# event clusters: (bp1 interval, bp2 interval) in alignment columns, from events_table_HM.tsv
events = {
 "A+B Henan (merged)":     ((1, 25), (578, 725)),
 "C EU344803/EU499311":    ((53, 71), (578, 681)),
 "EU344802":               ((1, 25), (578, 681)),
 "EU344805":               ((45, 53), (671, 681)),
 "OR134497":               ((1101, 1116), (1640, 1714)),
 "OR134498":               ((870, 973), (1646, 1655)),
 "OK094644":               ((197, 221), (581, 671)),
 "KP943594":               ((179, 191), (1116, 1129)),
}
sig = {  # one entry per breakpoint signature (group): groups 1,2(a),2(b),3,4; group 2 has two different genomes
 "1": ((13, 626)), "2a": ((1108, 1677)), "2b": ((921, 1650)), "3": ((209, 626)), "4": ((185, 1122))}
ev = np.array([[mid(*a), mid(*b)] for a, b in events.values()])
sg = np.array(list(sig.values()), dtype=float)

def maxwin(pos, win):
    pos = np.sort(pos % L)
    ext = np.concatenate([pos, pos + L])
    n = len(pos)
    best = 0
    j = 0
    for i in range(n):
        while ext[j] - ext[i] > win: j += 1
        # count of points in [ext[i], ext[i]+win]
    # vectorised: for each start i, count = searchsorted(ext, ext[i]+win, 'right') - i
    cnt = np.searchsorted(ext, pos + win, side="right") - np.arange(n)
    return cnt.max()

def maxwin_batch(P, win):
    P = np.sort(P % L, axis=1); n = P.shape[1]
    ext = np.concatenate([P, P + L], axis=1)
    best = np.zeros(P.shape[0], dtype=int)
    for i in range(n):
        c = (ext[:, i:i + n] <= (ext[:, [i]] + win)).sum(1)
        best = np.maximum(best, c)
    return best

recs = list(SeqIO.parse(sys.argv[1], "fasta"))
M = np.array([list(str(r.seq).upper()) for r in recs])
inf_cols = []
for j in range(M.shape[1]):
    col = M[:, j]; ok = np.isin(col, list("ACGT"))
    v, c = np.unique(col[ok], return_counts=True)
    if len(v) == 2 and c.min() >= 2: inf_cols.append(j + 1)
inf_cols = np.array(inf_cols); print(f"{len(inf_cols)} biallelic informative columns")

for label, pts in (("8 events, Henan A+B merged (16 breakpoints)", ev),):
    k = pts.shape[0]
    for win in (100, 150):
        obs = maxwin(pts.ravel(), win)
        u = rng.uniform(0, L, size=(NSIM, 2 * k))
        pu = (maxwin_batch(u, win) >= obs).mean()
        ic = rng.choice(inf_cols, size=(NSIM, 2 * k)).astype(float)
        pi = (maxwin_batch(ic, win) >= obs).mean()
        off = rng.uniform(0, L, size=(NSIM, k))
        rot = np.concatenate([pts[:, 0][None, :] + off, pts[:, 1][None, :] + off], axis=1)
        pr = (maxwin_batch(rot, win) >= obs).mean()
        print(f"{label}, window {win} nt: observed max = {obs}   p(uniform) = {pu:.4f}   p(informative sites) = {pi:.4f}   p(rotated arcs) = {pr:.4f}")
