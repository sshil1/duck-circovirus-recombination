import sys, json, math
from Bio import SeqIO
js = sys.argv[1] if len(sys.argv) > 1 else "gard_mpi2/gard_result.json"
raw = {r.id: str(r.seq).upper() for r in SeqIO.parse("setA_aln_raw.fasta", "fasta")}
tre = {r.id: str(r.seq).upper() for r in SeqIO.parse("tree_aln.fasta", "fasta")}
ids = list(raw); n = len(ids); L = len(raw[ids[0]])
common = [i for i in tre if i in raw][:30]
cols = None
for name, ch in (("gap", "-"), ("gap+N", "-N")):
    c = [j for j in range(L) if sum(raw[i][j] not in ch for i in ids) / n >= 0.95]
    ok = common and len(c) == len(next(iter(tre.values()))) and all("".join(raw[i][j] for j in c) == tre[i] for i in common)
    print(f"occupancy rule '{name}': {len(c)} cols  reproduces tree_aln: {bool(ok)}")
    if ok: cols = c; break
if cols is None: sys.exit("FAIL: could not reproduce tree_aln column set - stop and paste this output")
refid = next((i for i in raw if i.startswith("HM162346")), None)
if refid is None: refid = ids[0]; print("WARNING: HM162346 not in Set A; using", refid, "- ORF coordinates will not match")
ref = raw[refid]; coord = []; k = 0
for ch in ref:
    if ch != "-": k += 1
    coord.append(k)
d = json.load(open(js))
print("file:", js, " sequences:", d["input"]["number of sequences"], " sites:", d["input"]["number of sites"])
best = max(d["improvements"], key=int)
print("best model key:", best, " deltaAICc:", round(d["improvements"][best]["deltaAICc"], 1))
segs = [d["breakpointData"][k]["bps"][0] for k in sorted(d["breakpointData"], key=int)]
print("GARD segments (trimmed cols, 1-based):", segs)
REP, CAP = (48, 927), (1156, 1930)
rows = []
for a, b in segs:
    rows.append((a, b, coord[cols[a-1]], coord[cols[b-1]]))
print("\nSEGMENTS in genome coordinates (", refid, ")")
for i, (a, b, x, y) in enumerate(rows, 1): print(f"  seg{i}: trimmed {a}-{b}  genome nt {x}-{y}")
print("\nBREAKPOINTS")
out = ["bp\tlast_nt_seg_left\tfirst_nt_seg_right\tin_ORF"]
for i in range(len(rows) - 1):
    c = rows[i][1]; X = coord[cols[c-1]]; Y = coord[cols[c]]
    where = []
    if REP[0] <= X <= REP[1]:
        p = X - REP[0] + 1; where.append(f"Rep CDS nt {p} (codon {math.ceil(p/3)}, pos {(p-1)%3+1})")
    if CAP[0] <= X <= CAP[1]:
        p = CAP[1] - X + 1; where.append(f"Cap CDS nt {p} on minus strand (codon {math.ceil(p/3)}, pos {(p-1)%3+1})")
    w = "; ".join(where) or "intergenic/non-ORF"
    print(f"  bp{i+1}: between genome nt {X} and {Y}  -> {w}")
    out.append(f"{i+1}\t{X}\t{Y}\t{w}")
open("gard_breakpoints.tsv", "w").write("\n".join(out) + "\n")
print("\nnote: genome is circular, the alignment origin is an artificial end; wrote gard_breakpoints.tsv")
