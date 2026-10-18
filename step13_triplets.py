import numpy as np
from Bio import SeqIO
S = {r.id.split("|")[0]: str(r.seq).upper() for r in SeqIO.parse("tree_aln.fasta", "fasta")}
ref = S["HM162346.1"]; coord = []; k = 0
for ch in ref:
    if ch != "-": k += 1
    coord.append(k)
T = [("EU344802.1","GQ423742.1","DQ166837.1"),("EU344805.1","GQ423742.1","DQ166837.1"),("OK094644.1","GQ423742.1","FJ554673.1"),
     ("MN928811.1","GQ423745.1","MN928809.1"),("MN928792.1","GQ423745.1","MN928809.1"),("MN928806.1","GQ423745.1","MN928805.1"),
     ("MN928803.1","GQ423745.1","MN928807.1"),("MN928804.1","GQ423745.1","MT646346.1"),("MN928802.1","GQ423745.1","MN928807.1"),
     ("EU499311.1","GQ423743.1","GQ423742.1"),("EU344803.1","GQ423743.1","GQ423742.1"),
     ("OR134497.1","OR134499.1","OR134469.1"),("OR134498.1","OR134499.1","MN928808.1"),("KP943594.1","PX471097.1","KC851822.1")]
def maxchi(x, m=5):
    n = len(x); a = np.cumsum(x); tot = a[-1]; best, bk = 0.0, None
    for k in range(m, n - m + 1):
        l1 = a[k-1]; l0 = k - l1; r1 = tot - l1; r0 = n - k - r1
        den = (l1 + l0) * (r1 + r0) * (l1 + r1) * (l0 + r0)
        if den == 0: continue
        c = n * (l1 * r0 - l0 * r1) ** 2 / den
        if c > best: best, bk = c, k
    return best, bk
rng = np.random.default_rng(1)
print("child | P1 (left-ish) / P2 | informative sites | best cut (genome nt) | child=P1 left% / right% | maxchi | perm p (2000)")
for c, p1, p2 in T:
    C, A, B = S[c], S[p1], S[p2]; cols = []; x = []
    for j in range(len(C)):
        if A[j] in "ACGT" and B[j] in "ACGT" and A[j] != B[j] and C[j] in (A[j], B[j]):
            cols.append(j); x.append(1 if C[j] == A[j] else 0)
    x = np.array(x)
    if len(x) < 15: print(c, p1, p2, "too few informative sites:", len(x)); continue
    obs, bk = maxchi(x)
    cnt = sum(maxchi(rng.permutation(x))[0] >= obs for _ in range(2000))
    p = (1 + cnt) / 2001
    print(f"{c} | {p1} / {p2} | {len(x)} | {coord[cols[bk]]} | {100*x[:bk].mean():.0f}% / {100*x[bk:].mean():.0f}% | {obs:.1f} | {p:.4f}")
