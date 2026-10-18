#!/usr/bin/env bash
set -euo pipefail
JS=${1:-gard_B/gard_result.json}; OUT=${2:-segB}; T=${3:-12}
set +u; source /opt/conda/etc/profile.d/conda.sh; conda activate cavphylo; set -u
cd ~/dcv; mkdir -p $OUT
python3 - "$JS" "$OUT" <<'PY'
import sys, json
from Bio import SeqIO
js, out = sys.argv[1:3]
d = json.load(open(js))
segs = [d["breakpointData"][k]["bps"][0] for k in sorted(d["breakpointData"], key=int)]
recs = list(SeqIO.parse("tree_aln.fasta", "fasta"))
with open(f"{out}/segments.tsv", "w") as f:
    for i, (a, b) in enumerate(segs, 1):
        sub = [r[a-1:b] for r in recs]
        for r in sub: r.description = ""
        SeqIO.write(sub, f"{out}/seg{i}.fasta", "fasta")
        var = sum(1 for j in range(b-a+1) if len({str(r.seq[j]).upper() for r in sub if str(r.seq[j]) not in "-N"}) > 1)
        f.write(f"seg{i}\t{a}\t{b}\t{b-a+1}\t{var}\n")
print("seg\tfrom\tto\tcols\tvariable_sites"); print(open(f"{out}/segments.tsv").read())
PY
for s in $(ls $OUT/seg*.fasta | sort -V); do
  p=${s%.fasta}; echo "== $p $(date +%T)"
  iqtree3 -s $s -m MFP -B 1000 -bnni -alrt 1000 -nm 3000 -T $T --prefix $p -redo -quiet
done
echo ALL_DONE $(date +%T)
