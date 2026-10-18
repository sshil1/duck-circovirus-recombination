import numpy as np
from Bio import Phylo, SeqIO
meta = {}
for r in SeqIO.parse("tree_aln.fasta", "fasta"):
    f = r.id.split("|"); meta[f[0]] = f
accs = sorted(meta); idx = {a: i for i, a in enumerate(accs)}; N = len(accs)
exec(open("_td.py").read())
Ds = [tipdist(f"segB/seg{i}_f.contree") for i in (1, 2, 4)]
rows = [l.rstrip("\n").split("\t") for l in open("recomb_candidates.tsv")][1:]
flag = [idx[r[0]] for r in rows if float(r[3]) > 0.5]; fs = set(flag)
others = np.array([i for i in range(N) if i not in fs])
pairs = [np.sort(Ds[k][np.ix_(others, others)][np.triu_indices(len(others), 1)]) for k in range(3)]
print(len(flag), "extreme genomes (score>0.5)")
print("pairwise tree distance among the rest (median / 90th pct) seg1,seg2,seg4:", [(round(float(np.median(p)), 3), round(float(np.percentile(p, 90)), 3)) for p in pairs])
for k, nm in enumerate(["seg1", "seg2", "seg4"]):
    nnf = [int(np.argmin(np.where(np.arange(N) == g, np.inf, Ds[k][g]))) in fs for g in flag]
    print(f"{nm}: {sum(nnf)} of {len(flag)} have a flagged genome as nearest relative")
def rk(k, a, b): return np.searchsorted(pairs[k], Ds[k][a, b]) / len(pairs[k])
print("\ngenome country year | parent in seg1, seg2, seg4 | d(genome,parent) | parent-separation percentile: seg1v4(in seg4 tree) seg2v4 seg1v2")
for g in flag:
    par = [int(others[np.argmin(Ds[k][g, others])]) for k in range(3)]
    dg = [round(float(Ds[k][g, par[k]]), 3) for k in range(3)]
    m = meta[accs[g]]
    print(accs[g], m[3] if len(m) > 3 else "", m[4] if len(m) > 4 else "", "|", ",".join(accs[p] for p in par), "|", dg, "|", round(rk(2, par[0], par[2]), 2), round(rk(2, par[1], par[2]), 2), round(rk(1, par[0], par[1]), 2))
