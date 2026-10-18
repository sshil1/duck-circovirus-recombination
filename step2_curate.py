#!/usr/bin/env python3
"""STEP 2 - parse, flag, deduplicate, tier. Nothing is silently dropped: every record stays in dcv_curated.tsv with reasons.
Outputs: dcv_curated.tsv, dcv_core.fasta (original orientation), dcv_curation_report.txt
"""
import csv, hashlib, re, sys
from collections import Counter
from Bio import SeqIO
LMIN, LMAX = 1700, 2100   # DuCV-3 is ~1.75 kb; DuCV-1/2 ~1.99 kb; 'near-complete' records ~1.88 kb are kept but flagged
FULL = 1960
def patent(a):
    c = a.split(".")[0]
    return bool(re.match(r"^[AEI]\d{5}$", c) or re.match(r"^(DI|DJ|DL|DM|HH|HI|HV|HW|JA|JB|JC|JD|JE|PF|PI|CQ|CS|FV|GM|GN|LP|LQ|MA|MP)\d{6}$", c))
DUCK = re.compile(r"duck|anas|mallard|pekin|muscovy|cairina|teal|aythya|melanitta|scoter|mareca|wigeon|pintail|spatula|shoveler|tadorna|shelduck|eider|somateria|goldeneye|bucephala|merganser|mergus|aix|netta|pochard|dendrocygna|anatidae", re.I)
GOOSE = re.compile(r"goose|geese|anser|cygnopsis", re.I)
def hostclass(h, desc):
    t = f"{h} {desc}"
    if DUCK.search(h or ""): return "duck"
    if GOOSE.search(h or ""): return "goose"
    if re.search(r"swan|cygnus|gull|heron|egret|bird|avian|chicken|gallus|pigeon|columba|vulture|eagle", h or "", re.I): return "other_avian"
    if h: return "other_nonduck"
    if DUCK.search(desc or ""): return "duck_by_title"
    return "unknown"
rows = []; seqs = {}
for rec in SeqIO.parse("dcv_raw.gb", "genbank"):
    s = str(rec.seq).upper(); f = next((x for x in rec.features if x.type == "source"), None)
    q = f.qualifiers if f else {}; g = lambda k: (q.get(k) or [""])[0]
    geo = g("country") or g("geo_loc_name")
    country = geo.split(":")[0].strip()
    year = re.search(r"(19|20)\d{2}", g("collection_date")); 
    refs = [r for r in rec.annotations.get("references", [])]
    rows.append(dict(accession=rec.id, length=len(s), organism=rec.annotations.get("organism", ""), strain=g("strain") or g("isolate"),
        host=g("host"), completeness="", country=country, region=geo.partition(":")[2].strip(), collection_year=year.group(0) if year else "",
        description=rec.description, hostclass=hostclass(g("host"), rec.description),
        n_amb=sum(c not in "ACGT" for c in s), first_title=(refs[0].title if refs else ""), authors=(refs[0].authors[:60] if refs else ""),
        journal=(refs[0].journal[:60] if refs else ""), pubmed=(refs[0].pubmed_id if refs else ""), md5=hashlib.md5(s.encode()).hexdigest()))
    seqs[rec.id] = s
seen = {}
for r in rows:
    ex = []
    if not (LMIN <= r["length"] <= LMAX): ex.append(f"length {r['length']} outside {LMIN}-{LMAX}")
    if patent(r["accession"]): ex.append("patent-division sequence")
    if r["n_amb"] > 10: ex.append(">10 ambiguous bases")
    if re.search(r"vaccine|attenuat|infectious clone|clone", r["description"] + r["strain"], re.I): ex.append("vaccine/clone flag - check manually")
    if r["hostclass"] in ("other_nonduck",): ex.append(f"non-waterfowl host ({r['host']})")
    if r["md5"] in seen: ex.append(f"identical to {seen[r['md5']]}")
    else: seen[r["md5"]] = r["accession"]
    r["completeness"] = "full" if r["length"] >= FULL else "short(<%d)" % FULL
    r["exclude_reasons"] = "; ".join(ex); r["usable_core"] = "" if ex else "yes"
    r["usable_meta"] = "yes" if (r["usable_core"] and r["country"] and r["collection_year"]) else ""
cols = list(rows[0].keys())
with open("dcv_curated.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
with open("dcv_core.fasta", "w") as f:
    for r in rows:
        if r["usable_core"]:
            f.write(f">{r['accession']}|{(r['strain'] or 'NA').replace(' ','_')[:30]}|{r['hostclass']}|{r['country'] or 'NA'}|{r['collection_year'] or 'NA'}\n{seqs[r['accession']]}\n")
core = [r for r in rows if r["usable_core"]]
# submission blocks = runs of consecutive accession numbers with the same prefix
def split(a): m = re.match(r"([A-Z]+)(\d+)", a); return m.group(1), int(m.group(2))
srt = sorted(core, key=lambda r: split(r["accession"])); blocks = []; cur = []
for r in srt:
    if cur and split(r["accession"])[0] == split(cur[-1]["accession"])[0] and split(r["accession"])[1] - split(cur[-1]["accession"])[1] <= 1: cur.append(r)
    else:
        if cur: blocks.append(cur)
        cur = [r]
if cur: blocks.append(cur)
bid = {}
for i, b in enumerate(sorted(blocks, key=len, reverse=True)):
    for r in b: bid[r["accession"]] = i
for r in rows: r["block"] = bid.get(r["accession"], "")
with open("dcv_curated.tsv", "w", newline="") as f:
    cols = list(rows[0].keys()); w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
L = []; P = L.append
P(f"fetched {len(rows)}   usable_core {len(core)}   usable_meta(country+year) {sum(1 for r in rows if r['usable_meta'])}")
P("\nexclusion reasons:"); t = Counter(re.sub(r"\d+", "N", x) for r in rows for x in r["exclude_reasons"].split("; ") if x)
[P(f"  {v:>4}  {k}") for k, v in t.most_common()]
P("\ncompleteness (core):"); [P(f"  {k:14}{v}") for k, v in Counter(r["completeness"] for r in core).most_common()]
P("\nhost class (core):"); [P(f"  {k:14}{v}") for k, v in Counter(r["hostclass"] for r in core).most_common()]
P("\norganism field (core):"); [P(f"  {v:>4}  {k}") for k, v in Counter(r["organism"] for r in core).most_common()]
P("\ncountry (core):"); [P(f"  {k or 'NA':16}{v}") for k, v in Counter(r["country"] for r in core).most_common()]
P("\nyear (core):"); [P(f"  {k or 'NA':6}{v}") for k, v in sorted(Counter(r["collection_year"] for r in core).items())]
bs = sorted((len(b) for b in blocks), reverse=True)
P(f"\nsubmission blocks: {len(blocks)}; top-5 hold {sum(bs[:5])} ({100*sum(bs[:5])//len(core)}%), top-12 hold {sum(bs[:12])} ({100*sum(bs[:12])//len(core)}%)")
P("largest blocks: " + ", ".join(f"{b[0]['accession']}-{b[-1]['accession']} n={len(b)} {b[0]['country']}" for b in sorted(blocks, key=len, reverse=True)[:8]))
P(f"\nlength distribution (core): {sorted(Counter(r['length'] for r in core).items())}")
P(f"no PubMed link: {sum(1 for r in core if not r['pubmed'])}/{len(core)}")
P("\nlength-excluded records outside window (check these are not real genomes):")
for r in rows:
    if r["exclude_reasons"].startswith("length"): P(f"  {r['accession']} {r['length']} {r['description'][:70]}")
open("dcv_curation_report.txt", "w").write("\n".join(L)); print("\n".join(L))
