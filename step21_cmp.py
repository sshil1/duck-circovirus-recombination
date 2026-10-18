from Bio import Phylo
focal="MN928811.1 MN928792.1 MN928806.1 MN928803.1 MN928804.1 MN928802.1 EU344802.1".split()
T={k:Phylo.read(f"rec3_{k}.nwk","newick") for k in ("blockL","blockR")}
for k,t in T.items():
    tips=[x.name for x in t.get_terminals()]
    print(f"\n== {k}: {len(tips)} tips ==")
    for f in focal:
        d=sorted((t.distance(f,o),o) for o in tips if o not in focal)[:3]
        print(f, "->", ", ".join(f"{o} ({x:.3f})" for x,o in d))
