#!/usr/bin/env python3
"""Recombination-free negative control for the PHI test (Bruen et al. 2006) on the DuCV alignment.

Question: how often does a site-order-permutation PHI test call a recombination-FREE alignment
"significant" when the sequences evolve on the real tree under the fitted model and real-looking
rate variation?  The positive controls only show power; this measures the false-positive rate.

Two recombination-free simulation models, both on the real ML tree (dcv_iqtree.treefile), with
the TIM3+F+R4 parameters fitted by IQ-TREE (GTR parameterisation, identical matrix):
  A  "iid"     each site draws its rate category independently (the assumption behind PHI's null)
  B  "clustered" each site keeps the rank of its observed variability in the REAL alignment, so
               fast and slow sites sit at the same positions as in the real data, and the marginal
               rate distribution is still the fitted R4 (spatial rate structure preserved)
Model B is the harder test: if spatial rate structure alone makes PHI significant, the real
p = 0.001 cannot be trusted.

Usage:  python3 phi_negctl.py tree_aln.fasta dcv_iqtree.treefile OUTPREFIX [NREP=100] [NPERM=199] [NPROC=8]
Outputs: OUTPREFIX_summary.txt, OUTPREFIX_reps.tsv
"""
import sys, os, numpy as np, io
from multiprocessing import Pool
from Bio import SeqIO, Phylo

ALN, TREE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NREP = int(sys.argv[4]) if len(sys.argv) > 4 else 100
NPERM = int(sys.argv[5]) if len(sys.argv) > 5 else 199
NPROC = int(sys.argv[6]) if len(sys.argv) > 6 else 8
W = 100

# ---- fitted TIM3+F+R4 (from dcv_iqtree.iqtree) ----
pi = np.array([0.222, 0.2438, 0.2653, 0.2688])            # A C G T
R = {("A", "C"): 0.5891, ("A", "G"): 1.6439, ("A", "T"): 1.0,
     ("C", "G"): 0.5891, ("C", "T"): 2.8984, ("G", "T"): 1.0}
cat_w = np.array([0.7125, 0.1807, 0.09103, 0.01578])
cat_r = np.array([0.1651, 1.634, 4.407, 11.79])
B_ = "ACGT"
Q = np.zeros((4, 4))
for (a, b), v in R.items():
    i, j = B_.index(a), B_.index(b)
    Q[i, j] = v * pi[j]; Q[j, i] = v * pi[i]
np.fill_diagonal(Q, -Q.sum(1))
Q /= -(pi * np.diag(Q)).sum()                              # 1 expected substitution per unit time
# symmetric eigendecomposition of the reversible Q
sq = np.sqrt(pi)
Ssym = (sq[:, None] * Q) / sq[None, :]
lam, U = np.linalg.eigh(Ssym)
def P(t, r):
    return ((U * np.exp(lam * r * t)) @ U.T) * (1 / sq[:, None]) * sq[None, :]

# ---- real alignment and tree ----
recs = list(SeqIO.parse(ALN, "fasta"))
names = [r.id for r in recs]
M = np.array([list(str(r.seq).upper()) for r in recs])
N, L = M.shape
tree = Phylo.read(TREE, "newick")
tips = [t.name for t in tree.get_terminals()]
assert set(tips) == set(names), "tree tips and alignment names differ"
nodes = list(tree.find_clades(order="preorder"))
idx = {id(c): k for k, c in enumerate(nodes)}
parent = [-1] * len(nodes); blen = [0.0] * len(nodes)
for c in nodes:
    for ch in c.clades:
        parent[idx[id(ch)]] = idx[id(c)]; blen[idx[id(ch)]] = max(ch.branch_length or 0.0, 0.0)
tipnode = {c.name: idx[id(c)] for c in tree.get_terminals()}
order_rows = [tipnode[n] for n in names]

_map = {"A": 0, "C": 1, "G": 2, "T": 3}
REALNUM = np.vectorize(lambda x: _map.get(x, -1))(M).astype(np.int8)
MISSING = REALNUM < 0
BLSCALE = float(os.environ.get("BLSCALE", "1.0"))
blen = [b * BLSCALE for b in blen]
print(f"alignment {N} x {L}; tree nodes {len(nodes)}; branch-length scale {BLSCALE}", flush=True)

# per-site variability of the real alignment (fraction not equal to majority state)
fvar = np.zeros(L)
for j in range(L):
    col = M[:, j]; ok = np.isin(col, list("ACGT"))
    if ok.sum() == 0: continue
    v, c = np.unique(col[ok], return_counts=True)
    fvar[j] = 1 - c.max() / ok.sum()

# rate-category assignment by rank, preserving positional structure (model B)
def clustered_cats(rng):
    jitter = fvar + rng.random(L) * 1e-9                   # random tie-break among equal sites
    rank = np.argsort(np.argsort(jitter))                  # 0 = slowest
    q = (rank + 0.5) / L
    cum = np.cumsum(cat_w) / cat_w.sum()
    # slowest quantile -> slowest category (cat 0 is the slowest)
    return np.searchsorted(cum, q)

def simulate(cats, rng):
    """One recombination-free alignment on the real tree. Returns N x L int8 (0-3)."""
    # pre-computed cumulative transition matrices per branch for the 4 categories
    states = np.zeros((len(nodes), L), dtype=np.int8)
    states[0] = rng.choice(4, size=L, p=pi / pi.sum())
    for k in range(1, len(nodes)):
        p = parent[k]; t = blen[k]
        cum = np.stack([np.cumsum(P(t, r), axis=1) for r in cat_r])   # (4 cats, 4 from, 4 to)
        u = rng.random(L)
        c = cum[cats, states[p]]                                       # (L, 4)
        states[k] = np.minimum((u[:, None] > c).sum(1), 3)
    return states[order_rows]

def sites(Mnum):
    """Biallelic informative sites, same rule as phi_own.py. Mnum: N x L ints 0-3, -1 = missing."""
    cols = []
    for j in range(Mnum.shape[1]):
        col = Mnum[:, j]; ok = col >= 0
        v, c = np.unique(col[ok], return_counts=True)
        if len(v) == 2 and c.min() >= 2:
            cols.append(np.where(ok, (col != v[np.argmax(c)]).astype(np.int8), -1))
    return np.array(cols).T

def incompat(Bm):
    S = Bm.shape[1]; I = np.zeros((S, S), np.int8); valid = Bm >= 0
    for d in range(1, S):
        a, b, v = Bm[:, :-d], Bm[:, d:], valid[:, :-d] & valid[:, d:]
        inc = (((a == 0) & (b == 0) & v).any(0) & ((a == 0) & (b == 1) & v).any(0) &
               ((a == 1) & (b == 0) & v).any(0) & ((a == 1) & (b == 1) & v).any(0)).astype(np.int8)
        k = np.arange(S - d); I[k, k + d] = inc; I[k + d, k] = inc
    return I

def phi_test(Mnum, rng):
    Bm = sites(Mnum); S = Bm.shape[1]
    w = W if S >= 2 * W else max(10, S // 5)
    I = incompat(Bm)
    ii, jj = np.triu_indices(S, 1); near = (jj - ii) <= w
    a, b = ii[near], jj[near]
    f = lambda o: I[o[a], o[b]].mean()
    obs = f(np.arange(S))
    null = np.array([f(rng.permutation(S)) for _ in range(NPERM)])
    p = (np.sum(null <= obs) + 1) / (NPERM + 1)
    z = (obs - null.mean()) / null.std(ddof=1)
    return S, obs, null.mean(), z, p

def one(args):
    model, rep = args
    rng = np.random.default_rng(1000 * (1 if model == "iid" else 2) + rep)
    cats = (rng.choice(4, size=L, p=cat_w / cat_w.sum()) if model == "iid"
            else clustered_cats(rng))
    X = simulate(cats, rng)
    X = np.where(MISSING, -1, X)                           # same gaps/Ns as the real alignment
    const = float(np.mean([len(np.unique(X[X[:, j] >= 0, j])) <= 1 for j in range(L)]))
    S, obs, nm, z, p = phi_test(X, rng)
    return model, rep, S, const, obs, nm, z, p

if __name__ == "__main__":
    # observed value on the real data for the same code path
    Mnum = REALNUM
    rng = np.random.default_rng(1)
    realres = phi_test(Mnum, rng)
    print(f"REAL  S={realres[0]}  PHI={realres[1]:.5f}  null={realres[2]:.5f}  z={realres[3]:.2f}  p={realres[4]:.4f}", flush=True)
    realconst = float(np.mean([len(np.unique(Mnum[Mnum[:, j] >= 0, j])) <= 1 for j in range(L)]))

    jobs = [(m, r) for m in ("iid", "clustered") for r in range(NREP)]
    res = []
    with Pool(NPROC) as pool:
        for k, r in enumerate(pool.imap_unordered(one, jobs), 1):
            res.append(r)
            if k % 10 == 0: print(f"  {k}/{len(jobs)} done", flush=True)
    with open(OUT + "_reps.tsv", "w") as fh:
        fh.write("model\trep\tS\tconst_frac\tPHI\tnull_mean\tz\tp\n")
        for r in sorted(res):
            fh.write("\t".join(str(x) for x in r) + "\n")
    lines = [f"REAL alignment: S={realres[0]}  constant-site fraction={realconst:.3f}  z={realres[3]:.2f}  p={realres[4]:.4f}",
             f"{NREP} recombination-free replicates per model, {NPERM} permutations each (min p = {1/(NPERM+1):.4f})", ""]
    for m in ("iid", "clustered"):
        a = np.array([r[2:] for r in res if r[0] == m], dtype=float)
        S, cf, obs, nm, z, p = a.T
        lines += [f"model {m}:",
                  f"  informative sites S: mean {S.mean():.0f} (range {S.min():.0f}-{S.max():.0f}); constant-site fraction {cf.mean():.3f}",
                  f"  z: mean {z.mean():.2f}, sd {z.std(ddof=1):.2f}, 5th pct {np.percentile(z,5):.2f}, min {z.min():.2f}",
                  f"  false-positive rate at p<0.05: {np.mean(p<0.05):.3f}   at p<0.01: {np.mean(p<=0.01):.3f}",
                  f"  replicates with z <= real z ({realres[3]:.2f}): {int(np.sum(z <= realres[3]))} of {len(z)}", ""]
    open(OUT + "_summary.txt", "w").write("\n".join(lines))
    print("\n".join(lines))
