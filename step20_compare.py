import csv
rec=list(csv.DictReader(open("seq3/3s.rec.csv")))
rc=[l.rstrip("\n").split("\t") for l in open("recomb_candidates.tsv")]
print("recomb_candidates.tsv header:", rc[0])
flag=[r[0] for r in rc[1:] if float(r[3])>0.5]
byc={r["C_ACCNUM"]:r for r in rec}
print("\n3SEQ best triplet for the %d flagged genomes:"%len(flag))
print("child\tparentP\tparentQ\tmin_rec_len\tbreakpoints\tP_flagged?\tQ_flagged?")
fs=set(flag); both=0
for c in flag:
    r=byc.get(c)
    if not r: print(c,"not in 3SEQ output"); continue
    pf=r["P_ACCNUM"] in fs; qf=r["Q_ACCNUM"] in fs
    both+= (pf or qf)
    print(c,r["P_ACCNUM"],r["Q_ACCNUM"],r["min_rec_length"],r["breakpoints"],pf,qf,sep="\t")
print("\nflagged genomes with a flagged parent:",both,"of",len(flag))
others=[r for r in rec if r["C_ACCNUM"] not in fs]
print("non-flagged children with a triplet:",len(others),"| median min_rec_length:",sorted(int(r["min_rec_length"]) for r in others)[len(others)//2])
