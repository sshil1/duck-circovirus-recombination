#!/usr/bin/env python3
"""STEP 1 - fetch every plausible duck circovirus (DuCV) record from NCBI as GenBank.
Three queries are unioned because titles/taxonomy are inconsistent (DuCV-3 may sit under 'Circovirus sp.'):
  A  txid324685 subtree, 1500-2500 bp
  B  title contains duck circovirus / DuCV, 1500-2500 bp (any taxid)
  C  title contains 'duck circovirus' / 'DuCV', 1500-2500 bp, organism 'Circovirus' subtree (catches mis-taxonomised)
Outputs: dcv_raw.gb, uids.json.   Usage: python3 step1_fetch.py you@example.org [NCBI_API_KEY]
"""
import json, sys, time, urllib.parse, urllib.request
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = sys.argv[1] if len(sys.argv) > 1 else "none@example.org"
KEY = sys.argv[2] if len(sys.argv) > 2 else None
def eget(tool, params, tries=6):
    params = dict(params, tool="dcv_kit", email=EMAIL)
    if KEY: params["api_key"] = KEY
    url = EU + tool + "?" + urllib.parse.urlencode(params); last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read()
        except Exception as e: last = e; time.sleep(3 * (i + 1))
    raise RuntimeError(f"{tool}: {last}")
def esearch(term):
    j = json.loads(eget("esearch.fcgi", {"db": "nuccore", "term": term, "retmode": "json", "retmax": 10000}))
    time.sleep(0.4); return j["esearchresult"]["idlist"]
Q = {"A": 'txid324685[Organism:exp] AND 1500:2500[SLEN]',
     "B": '("duck circovirus"[Title] OR DuCV[Title]) AND 1500:2500[SLEN]',
     "C": '(duck[Title] AND circovirus[Title]) AND txid39725[Organism:exp] AND 1500:2500[SLEN]'}
res = {k: esearch(v) for k, v in Q.items()}
for k in res: print(f"query {k}: {len(res[k])}")
union = sorted(set().union(*res.values()), key=int)
print(f"union: {len(union)}   (A-only {len(set(res['A'])-set(res['B'])-set(res['C']))}, outside A: {len(set(union)-set(res['A']))})")
json.dump({"queries": Q, **res, "union": union}, open("uids.json", "w"))
with open("dcv_raw.gb", "wb") as out:
    for i in range(0, len(union), 150):
        out.write(eget("efetch.fcgi", {"db": "nuccore", "id": ",".join(union[i:i+150]), "rettype": "gb", "retmode": "text"}))
        print(f"  fetched {min(i+150,len(union))}/{len(union)}", flush=True); time.sleep(0.5)
print("done -> dcv_raw.gb")
