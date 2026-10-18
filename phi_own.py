#!/usr/bin/env python3
"""Independent PHI (pairwise homoplasy index) test - Bruen, Philippe & Bryant (2006).

Why not PhiPack: its -o extras (NSS, Max Chi^2) are all-pairs scans that do not finish on a
551-taxon alignment. The PHI statistic itself is cheap if you precompute the site-by-site
incompatibility matrix ONCE and then permute site ORDER, which is exactly what the null
hypothesis does - incompatibility values are fixed, only window membership changes.

Incompatibility = four-gamete test on biallelic informative sites. Under no recombination,
incompatible site pairs (homoplasies) are independent of how far apart the sites are. Under
recombination, nearby sites are MORE compatible than distant ones, so the windowed mean
incompatibility falls below the permuted null.

Null: 1000 random permutations of site order. One-sided (recombination lowers PHI).
"""
import sys, numpy as np
from Bio import SeqIO

W = 100          # window, in units of informative sites (PhiPack default)
NPERM = 1000
rng = np.random.default_rng(20261018)

recs = list(SeqIO.parse(sys.argv[1] if len(sys.argv) > 1 else "tree_aln.fasta", "fasta"))
M = np.array([list(str(r.seq).upper()) for r in recs])
N, L = M.shape
print(f"alignment: {N} taxa x {L} columns")

# biallelic informative sites: exactly 2 ACGT states, each in >=2 taxa
cols, B = [], []
for j in range(L):
    col = M[:, j]
    ok = np.isin(col, list("ACGT"))
    vals, cnts = np.unique(col[ok], return_counts=True)
    if len(vals) == 2 and cnts.min() >= 2:
        major = vals[np.argmax(cnts)]
        b = np.where(ok, (col != major).astype(np.int8), -1)   # -1 = missing
        cols.append(j); B.append(b)
B = np.array(B).T                                               # N x S
S = B.shape[1]
print(f"biallelic informative sites: {S}")
if S < 2 * W:
    W = max(10, S // 5); print(f"  (window reduced to {W})")

# incompatibility matrix via the four-gamete test, computed by offset
valid = B >= 0
I = np.zeros((S, S), dtype=np.int8)
for d in range(1, S):
    a, b = B[:, :-d], B[:, d:]
    v = valid[:, :-d] & valid[:, d:]
    g00 = ((a == 0) & (b == 0) & v).any(axis=0)
    g01 = ((a == 0) & (b == 1) & v).any(axis=0)
    g10 = ((a == 1) & (b == 0) & v).any(axis=0)
    g11 = ((a == 1) & (b == 1) & v).any(axis=0)
    inc = (g00 & g01 & g10 & g11).astype(np.int8)               # all four gametes -> incompatible
    idx = np.arange(S - d)
    I[idx, idx + d] = inc; I[idx + d, idx] = inc
print(f"incompatible site pairs: {int(I.sum()//2)} of {S*(S-1)//2} "
      f"({100*I.sum()/ (S*(S-1)):.1f}%)")

# windowed mean incompatibility under a given site ordering
ii, jj = np.triu_indices(S, 1)
near = (jj - ii) <= W
ii_n, jj_n = ii[near], jj[near]
def phi(order): return I[order[ii_n], order[jj_n]].mean()

obs = phi(np.arange(S))
null = np.array([phi(rng.permutation(S)) for _ in range(NPERM)])
p = (np.sum(null <= obs) + 1) / (NPERM + 1)                      # one-sided
z = (obs - null.mean()) / null.std(ddof=1)

print(f"\nPHI observed (mean incompatibility, window={W}): {obs:.5f}")
print(f"permuted null: mean {null.mean():.5f}  sd {null.std(ddof=1):.5f}  "
      f"range {null.min():.5f}-{null.max():.5f}")
print(f"z = {z:.2f}")
print(f"one-sided permutation p = {p:.4f}  ({NPERM} permutations)")
print("\nVERDICT:", "SIGNIFICANT RECOMBINATION (p < 0.05)" if p < 0.05 else
      "no significant recombination signal (p >= 0.05)")
