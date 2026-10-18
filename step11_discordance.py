import sys, itertools, numpy as np
from Bio import Phylo, SeqIO
segs = [f"segB/seg{i}_f.contree" for i in range(1, 7)]
meta = {}
for r in SeqIO.parse("tree_aln.fasta", "fasta"):
    f = r.id.split("|"); meta[f[0]] = f
accs = sorted(meta); idx = {a: i for i, a in enumerate(accs)}; N = len(accs)
def tipdist(path):
    t = Phylo.read(path, "newick"); tips = t.get_terminals()
    parent = {}; depth = {t.root: 0.0}
    for cl in t.find_clades(order="level"):
        for ch in cl.clades:
            parent[ch] = cl; depth[ch] = depth[cl] + (ch.branch_length or 0.0)
    anc = {}
    for tip in tips:
        a, n = [], tip
        while True:
            a.append(n)
            if n in parent: n = parent[n]
            else: break
        anc[tip] = a
    ids = [idx[tp.name.split("|")[0]] for tp in tips]
    D = np.zeros((N, N)); sets = [set(anc[tp]) for tp in tips]
    for x in range(len(tips)):
        for y in range(x + 1, len(tips)):
            for n in anc[tips[y]]:
                if n in sets[x]: l = depth[n]; break
            d = depth[tips[x]] + depth[tips[y]] - 2 * l
            D[ids[x], ids[y]] = D[ids[y], ids[x]] = d
    return D
Ds = []
for s in segs:
    Ds.append(tipdist(s)); print("tree", s, "done", flush=True)
def score(pairs):
    out = np.zeros(N)
    for g in range(N):
        m = np.arange(N) != g; cs = []
        for i, j in pairs:
            a, b = Ds[i][g, m], Ds[j][g, m]
            cs.append(np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else np.nan)
        out[g] = 1 - np.nanmean(cs)
    return out
main = score(list(itertools.combinations([0, 1, 3], 2)))
allp = score(list(itertools.combinations(range(6), 2)))
med = np.median(main); mad = 1.4826 * np.median(np.abs(main - med)); thr = med + 3 * mad
nn = [[accs[int(np.argmin(np.where(np.arange(N) == g, np.inf, Ds[s][g])))] for s in range(6)] for g in range(N)]
order = np.argsort(-main)
with open("recomb_candidates.tsv", "w") as f:
    f.write("acc\tcountry\tyear\tscore_main(seg1,2,4)\tscore_all\tflag\tNN_seg1..6\n")
    for g in order:
        m = meta[accs[g]]
        f.write(f"{accs[g]}\t{m[3] if len(m)>3 else ''}\t{m[4] if len(m)>4 else ''}\t{main[g]:.4f}\t{allp[g]:.4f}\t{'FLAG' if main[g] > thr else ''}\t{','.join(nn[g])}\n")
print(f"median {med:.4f}  MAD-sigma {mad:.4f}  threshold {thr:.4f}  flagged {int((main > thr).sum())} of {N}")
print("top 20:")
for g in order[:20]:
    m = meta[accs[g]]; print(" ", accs[g], m[3] if len(m) > 3 else "", m[4] if len(m) > 4 else "", round(main[g], 3), "NN:", ",".join(nn[g][:4]))
print("wrote recomb_candidates.tsv")
