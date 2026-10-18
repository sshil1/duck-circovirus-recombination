#!/usr/bin/env python3
"""Positive control + distance-decay profile for the PHI test.

A negative PHI is only meaningful if the test can detect recombination that IS present.
Control: take the real alignment and manufacture mosaics - for a fraction of taxa, replace a
contiguous block with the same block from a distantly related taxon. That is exactly the
signal PHI is built to find. If PHI stays non-significant on manufactured recombinants, the
test lacks power here and the real negative means nothing.
"""
import numpy as np
from Bio import SeqIO
rng = np.random.default_rng(7)
W, NPERM = 100, 1000

def sites(M):
    N, L = M.shape; B = []
    for j in range(L):
        col = M[:, j]; ok = np.isin(col, list("ACGT"))
        v, c = np.unique(col[ok], return_counts=True)
        if len(v) == 2 and c.min() >= 2:
            B.append(np.where(ok, (col != v[np.argmax(c)]).astype(np.int8), -1))
    return np.array(B).T

def incmat(B):
    S = B.shape[1]; valid = B >= 0; I = np.zeros((S, S), np.int8)
    for d in range(1, S):
        a, b, v = B[:, :-d], B[:, d:], valid[:, :-d] & valid[:, d:]
        inc = (((a==0)&(b==0)&v).any(0) & ((a==0)&(b==1)&v).any(0) &
               ((a==1)&(b==0)&v).any(0) & ((a==1)&(b==1)&v).any(0)).astype(np.int8)
        idx = np.arange(S-d); I[idx, idx+d] = inc; I[idx+d, idx] = inc
    return I

def phitest(M, label):
    B = sites(M); S = B.shape[1]; I = incmat(B)
    ii, jj = np.triu_indices(S, 1); near = (jj-ii) <= W
    f = lambda o: I[o[ii[near]], o[jj[near]]].mean()
    obs = f(np.arange(S))
    null = np.array([f(rng.permutation(S)) for _ in range(NPERM)])
    p = (np.sum(null <= obs)+1)/(NPERM+1)
    z = (obs-null.mean())/null.std(ddof=1)
    print(f"  {label:34} S={S:4}  PHI={obs:.5f}  null={null.mean():.5f}  z={z:+6.2f}  p={p:.4f}  "
          f"{'SIGNIF' if p<0.05 else 'ns'}")
    # distance decay
    if "real" in label:
        print("\n  incompatibility by site separation (real data):")
        for lo, hi in ((1,20),(21,50),(51,100),(101,170),(171,S-1)):
            m = (jj-ii >= lo) & (jj-ii <= hi)
            if m.sum(): print(f"    sep {lo:>4}-{hi:<4} sites: {I[ii[m], jj[m]].mean():.4f}  (n={m.sum()})")
        print()
    return p

M = np.array([list(str(r.seq).upper()) for r in SeqIO.parse("tree_aln.fasta", "fasta")])
N, L = M.shape
print(f"alignment {N} x {L}\n")
phitest(M, "real alignment")

print("positive controls - manufactured recombinants:")
for frac in (0.05, 0.10, 0.25, 0.50):
    R = M.copy()
    n = int(N*frac)
    recip = rng.choice(N, n, replace=False)
    for r in recip:
        donor = rng.integers(0, N)
        s = rng.integers(int(L*0.2), int(L*0.8))        # one breakpoint, block to the 3' end
        R[r, s:] = M[donor, s:]
    phitest(R, f"{int(frac*100)}% taxa recombinant")
