import json, os
os.chdir(os.path.expanduser("~/dcv/sel"))
OFF = {"repA": 0, "repB": 210, "capA": 0, "capB": 78}
def refmap(P):
    seq, cur = {}, None
    for l in open(f"{P}.fasta"):
        l = l.strip()
        if l.startswith(">"): cur = l[1:]; seq[cur] = []
        elif cur: seq[cur].append(l)
    r = "".join(seq.get("HM162346_1", []))
    if not r: return None
    out, c = [], 0
    for i in range(0, len(r), 3):
        if r[i:i+3] == "---": out.append(None)
        else: c += 1; out.append(c + OFF[P])
    return out
def load(P, m):
    f = f"hyphy/{P}.{m}.json"
    try: return json.load(open(f))
    except Exception: return None
def hs(d): return [h[0] for h in d["MLE"]["headers"]]
def rows(d): return d["MLE"]["content"]["0"]
rowsout = []
for P in ["repA", "repB", "capA", "capB"]:
    rm = refmap(P); line = {"part": P}
    d = load(P, "busted")
    if d:
        tr = d['test results']; pv = tr.get('p', tr.get('p-value'))
        line['BUSTED'] = 'LRT=%s p=%s' % (round(tr.get('LRT', float('nan')), 2), pv)
        if pv is None: print('BUSTED keys:', list(tr))
    sets = {}
    d = load(P, "fel")
    if d:
        h = hs(d); ip = h.index("p-value"); R = rows(d)
        pos = {i for i, r in enumerate(R) if r[ip] <= 0.1 and r[1] > r[0]}
        neg = {i for i, r in enumerate(R) if r[ip] <= 0.1 and r[1] < r[0]}
        sets["FEL+"] = pos; line["FEL_pos/neg"] = f"{len(pos)}/{len(neg)}"
    d = load(P, "meme")
    if d:
        h = hs(d); ip = h.index("p-value"); R = rows(d)
        pos = {i for i, r in enumerate(R) if r[ip] <= 0.1}
        sets["MEME"] = pos; line["MEME_sites"] = str(len(pos))
    d = load(P, "fubar")
    if d:
        h = hs(d); ia = [i for i, x in enumerate(h) if "Prob[alpha<beta]" in x][0]; ib = [i for i, x in enumerate(h) if "Prob[alpha>beta]" in x][0]; R = rows(d)
        pos = {i for i, r in enumerate(R) if r[ia] >= 0.9}; neg = {i for i, r in enumerate(R) if r[ib] >= 0.9}
        sets["FUBAR+"] = pos; line["FUBAR_pos/neg"] = f"{len(pos)}/{len(neg)}"
    print(P, line)
    cnt = {}
    for k, s in sets.items():
        for i in s: cnt[i] = cnt.get(i, []) + [k]
    multi = sorted((i, v) for i, v in cnt.items() if len(v) >= 2)
    if multi:
        print("  sites supported by >=2 positive-selection methods (alignment codon -> reference CDS codon):")
        for i, v in multi: print("   ", i + 1, "->", (rm[i] if rm and rm[i] else "insertion/NA"), v)
    else: print("  no site is supported by two or more positive-selection methods")
