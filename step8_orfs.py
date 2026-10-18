#!/usr/bin/env python3
"""STEP 8 - extract Rep (ORF1) and Cap (ORF2) coding sequences from every genome, QC them, and build codon-aware alignments.
Independent of GARD. Uses the step-4 nucleotide alignment (dcv_aln.fasta) and the reference's ORF coordinates (found on HM162346 by ORF scan):
  Rep : forward strand, genome positions [48, 927)    -> 292 aa + stop
  Cap : reverse strand, genome positions [1156, 1930) -> 257 aa + stop   (reverse-complemented)
ORF3 is NOT handled here (it overlaps other frames; needs its own treatment).
Per genome: take the alignment columns spanning the reference ORF, remove gaps, translate, QC (in-frame, starts ATG, ends with stop, no internal stop).
QC-pass CDSs -> protein MAFFT alignment -> back-translated to a codon alignment (frame-safe, like pal2nal).
Outputs: orf_qc.tsv, rep_codon_aln.fasta, cap_codon_aln.fasta, rep_aa_aln.fasta, cap_aa_aln.fasta, orf_summary.txt
Usage: python3 step8_orfs.py [threads]
"""
import csv, subprocess, sys
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
T = sys.argv[1] if len(sys.argv) > 1 else "4"
REF = "HM162346.1"
ORFS = {"rep": (48, 927, "+"), "cap": (1156, 1930, "-")}
aln = [(r.id, str(r.seq).upper()) for r in SeqIO.parse("dcv_aln.fasta", "fasta")]
ref_id, ref = next(x for x in aln if x[0].startswith(REF))
gpos2col = [c for c, ch in enumerate(ref) if ch != "-"]            # genome position -> alignment column
ref_g = ref.replace("-", ""); L = len(ref_g)
assert L == 1993, L
for name, (a, b, st) in ORFS.items():                                # sanity: the reference ORF must translate cleanly
    s = ref_g[a:b]; s = s if st == "+" else str(Seq(s).reverse_complement())
    p = str(Seq(s).translate()); assert p.startswith("M") and p.endswith("*") and "*" not in p[:-1], (name, p[:5], p[-5:])
    print(f"reference {name}: {len(p)-1} aa, OK")
codes = {}; qc = []
for rid, row in aln:
    acc = rid.split("|")[0]
    for name, (a, b, st) in ORFS.items():
        c0, c1 = gpos2col[a], gpos2col[b - 1]
        seg = row[c0:c1 + 1].replace("-", "")
        seg = seg if st == "+" else str(Seq(seg).reverse_complement())
        status = []
        if len(seg) % 3: status.append("frameshift/len%3")
        else:
            p = str(Seq(seg).translate())
            if not p.startswith("M"): status.append("no_ATG")
            if not p.endswith("*"): status.append("no_stop")
            if "*" in p[:-1]: status.append("internal_stop")
        if any(ch not in "ACGT" for ch in seg): status.append("ambiguous")
        gaps_edge = (row[c0] == "-" or row[c1] == "-")
        if gaps_edge: status.append("incomplete_end")
        qc.append((acc, rid, name, len(seg), ";".join(status) or "PASS"))
        if not status: codes.setdefault(name, []).append((rid, seg[:-3]))      # drop stop codon
with open("orf_qc.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["accession", "id", "orf", "cds_len", "qc"]); w.writerows(qc)
out = []
for name in ORFS:
    recs = codes.get(name, []); fa = f"{name}_aa.fasta"
    with open(fa, "w") as o:
        for rid, cds in recs: o.write(f">{rid}\n{Seq(cds).translate()}\n")
    with open(f"{name}_aa_aln.fasta", "w") as o: subprocess.run(["mafft", "--auto", "--thread", T, "--quiet", fa], stdout=o, check=True)
    cds = dict(recs); n_ok = 0
    with open(f"{name}_codon_aln.fasta", "w") as o:
        for r in SeqIO.parse(f"{name}_aa_aln.fasta", "fasta"):
            nt = cds[r.id]; k = 0; line = []
            for ch in str(r.seq):
                if ch == "-": line.append("---")
                else: line.append(nt[3*k:3*k+3]); k += 1
            assert k * 3 == len(nt); o.write(f">{r.id}\n{''.join(line)}\n"); n_ok += 1
    t = Counter(q[4] for q in qc if q[2] == name)
    out.append(f"{name}: {n_ok} codon-aligned of {len(aln)}; QC: {dict(t.most_common())}; codon alignment length {len(line)} codons? -> {len(''.join(line))//3}")
open("orf_summary.txt", "w").write("\n".join(out)); print("\n".join(out))
