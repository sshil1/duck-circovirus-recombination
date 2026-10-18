from Bio import Phylo
R6="MN928811.1 MN928792.1 MN928806.1 MN928803.1 MN928804.1 MN928802.1".split()
A=["GQ423745.1","GQ423743.1"]                                   # GQ4237xx group (left-block neighbours)
D=["MN928807.1","MN928808.1","MN928809.1","MN928805.1"]         # MN9288xx donors (right-block neighbours)
for k in ("blockL","blockR"):
    t=Phylo.read(f"iq_{k}.treefile","newick"); t.root_at_midpoint()
    tips={x.name for x in t.get_terminals()}
    miss=[n for n in R6+A+D if n not in tips]
    print(f"\n== {k} ({len(tips)} tips) missing names: {miss or 'none'} ==")
    for label,grp,other,oname in [("6 recombinants + GQ4237xx",A,D,"MN9288xx donors"),
                                   ("6 recombinants + MN9288xx donors",D,A,"GQ4237xx")]:
        c=t.common_ancestor(R6+grp)
        names={x.name for x in c.get_terminals()}
        inc=[o for o in other if o in names]
        print(f"{label}: clade of {len(names)} tips | support={c.confidence} | {oname} inside clade: {inc or 'none'}")
