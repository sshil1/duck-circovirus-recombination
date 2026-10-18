from Bio import SeqIO
bad = {l.split("\t")[0].split(".")[0] for l in list(open("recomb_candidates.tsv"))[1:] if float(l.split("\t")[3]) > 0.5}
REF = "HM162346.1"
def load(fn): return {r.id.split("|")[0]: str(r.seq).upper() for r in SeqIO.parse(fn, "fasta")}
def split(fn, boundary_nt, names):
    S = load(fn); ref = S[REF]; cnt = 0; cut = None
    for j, ch in enumerate(ref):
        if ch != "-": cnt += 1
        if cnt == boundary_nt: cut = j + 1; break
    cut = ((cut + 2) // 3) * 3
    print(f"{fn}: {len(S)} seqs, length {len(ref)}; boundary after reference CDS nt {boundary_nt} -> alignment column {cut} (codon {cut//3})")
    for nm, (a, b) in zip(names, [(0, cut), (cut, len(ref))]):
        seqs = {k: v[a:b] for k, v in S.items() if k.split(".")[0] not in bad}
        seen = {}
        for k in sorted(seqs, key=lambda x: (x != REF, x)):
            seen.setdefault(seqs[k], k)
        keep = sorted(seen.values())
        stops = sum(1 for k in keep for i in range(0, len(seqs[k]) - 2, 3) if seqs[k][i:i+3] in ("TAA", "TAG", "TGA"))
        var = sum(1 for j in range(b - a) if len({seqs[k][j] for k in keep if seqs[k][j] in "ACGT"}) > 1)
        with open(f"{nm}.fasta", "w") as o:
            for k in keep: o.write(f">{k.replace('.', '_')}\n{seqs[k]}\n")
        print(f"  {nm}: {(b-a)//3} codons, {len(seqs)} seqs -> {len(keep)} unique, variable sites {var}, in-frame stop codons {stops}")
split("rep_codon_aln.fasta", 630, ["repA", "repB"])
split("cap_codon_aln.fasta", 234, ["capA", "capB"])
