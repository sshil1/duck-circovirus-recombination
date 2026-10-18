import collections, numpy as np
from Bio import Phylo, SeqIO
meta = {}
for r in SeqIO.parse("tree_aln.fasta", "fasta"):
    f = r.id.split("|"); meta[f[0]] = f
accs = sorted(meta); idx = {a: i for i, a in enumerate(accs)}; N = len(accs)
exec(open("_td.py").read())
refs = {"1a": "MF627690", "1b": "HM162351", "1c": "GU014543", "1d": "KR491946", "2a": "EU344805", "2b": "KP229377", "2c": "ON227555"}
b2f = {a.split(".")[0]: a for a in accs}
present = {k: idx[b2f[v]] for k, v in refs.items() if v in b2f}
print("reference genomes found in Set A:", {k: b2f[refs[k]] for k in present}, "| MISSING:", [k for k in refs if k not in present])
names = ["whole", "s1", "s2", "s3", "s4", "s5", "s6"]
files = ["dcv_iqtree.treefile"] + [f"segB/seg{i}_f.contree" for i in range(1, 7)]
keys = list(present); cols = [present[k] for k in keys]; labs = {}
for nm, f in zip(names, files):
    D = tipdist(f); sub = D[:, cols]
    labs[nm] = [keys[int(np.argmin(sub[g]))] for g in range(N)]
    print(nm, "done", flush=True)
disc = [len({labs["s1"][g], labs["s2"][g], labs["s4"][g]}) > 1 for g in range(N)]
with open("genotype_by_segment.tsv", "w") as o:
    o.write("acc\tcountry\tyear\t" + "\t".join(names) + "\tdiscordant_s1_s2_s4\n")
    for g, a in enumerate(accs):
        m = meta[a]
        o.write(f"{a}\t{m[3] if len(m)>3 else ''}\t{m[4] if len(m)>4 else ''}\t" + "\t".join(labs[n][g] for n in names) + f"\t{disc[g]}\n")
print("whole-genome nearest-reference genotype counts:", sorted(collections.Counter(labs["whole"]).items()))
print("genomes with different genotype among seg1, seg2, seg4:", sum(disc), "of", N)
print("(seg1, seg2, seg4) combinations among those:")
for k, v in collections.Counter((labs["s1"][g], labs["s2"][g], labs["s4"][g]) for g in range(N) if disc[g]).most_common(12): print("  ", k, v)
print("\nthe 14 extreme genomes from step 11:")
for l in list(open("recomb_candidates.tsv"))[1:15]:
    a = l.split("\t")[0]; g = idx[a]
    print(" ", a, "whole:", labs["whole"][g], "| s1,s2,s4:", labs["s1"][g], labs["s2"][g], labs["s4"][g])
print("wrote genotype_by_segment.tsv")
