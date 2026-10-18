import collections
from Bio.Seq import Seq
gt = {}
for l in list(open("genotype_by_segment.tsv"))[1:]:
    f = l.rstrip("\n").split("\t"); gt[f[0]] = f[3]
def load(P):
    d, cur = {}, None
    for l in open(f"sel/{P}.fasta"):
        l = l.strip()
        if l.startswith(">"): cur = l[1:]; d[cur] = []
        elif cur: d[cur].append(l)
    return {k: "".join(v) for k, v in d.items()}
def acc(k): i = k.rfind("_"); return k[:i] + "." + k[i+1:]
SITES = [("capA", 5, 5), ("capA", 12, 12), ("capA", 18, 18), ("capA", 23, 23), ("capA", 49, 49), ("capB", 28, 106), ("repA", 73, 71)]
for P, i, ref in SITES:
    S = load(P); aas = {}
    for k, s in S.items():
        c = s[3*(i-1):3*i]
        aas[k] = str(Seq(c).translate()) if set(c) <= set("ACGT") else "-"
    cnt = collections.Counter(aas.values())
    print(f"{P} codon {i} (reference codon {ref}): {len(S)} unique seqs; residues: " + ", ".join(f"{a}:{n}" for a, n in cnt.most_common()))
    for a, n in cnt.most_common(4):
        if a == "-": continue
        g = collections.Counter(gt.get(acc(k), "?") for k, v in aas.items() if v == a)
        print(f"    {a} ({n}): genotypes " + ", ".join(f"{x}:{y}" for x, y in g.most_common(5)))
