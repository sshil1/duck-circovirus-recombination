import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, Patch
plt.rcParams["hatch.linewidth"]=0.6
rows=list(csv.DictReader(open("events_table_HM.tsv"),delimiter="\t"))
GEN=1993
SHORT={"1":"1: shared rep (~578-680) /\norigin-junction signature","2":"2: intergenic to cap","3":"3: within rep","4":"4: rep to intergenic\n(poor fit)"}
iv=lambda s:tuple(int(x) for x in s.split("-"))
OR,BL,GR,PK="#E69F00","#0072B2","#009E73","#CC79A7"
y=1; ys=[]; prev=None; gy={}
for r in rows:
    g=r["group"].split()[0]
    if prev is not None and g!=prev: y+=0.7
    ys.append(y); gy.setdefault(g,[]).append(y); y+=1; prev=g
fig,ax=plt.subplots(figsize=(12.5,0.34*len(rows)+3.2))
for a,b in ((928,1155),(1931,GEN),(1,47)): ax.add_patch(Rectangle((a,-0.6),b-a,y+0.3,color="#eeeeee",zorder=0))
ax.add_patch(FancyArrow(48,0,879,0,width=0.5,head_width=0.7,head_length=40,length_includes_head=True,color=GR))
ax.add_patch(FancyArrow(1930,0,-774,0,width=0.5,head_width=0.7,head_length=40,length_includes_head=True,color=PK))
ax.text(480,0,"rep",ha="center",va="center",color="white",fontsize=9)
ax.text(1545,0,"cap",ha="center",va="center",color="black",fontsize=9)
ax.text(1043,0,"intergenic",ha="center",va="center",fontsize=7,color="#555555")
ax.text(1964,-0.75,"intergenic",ha="center",va="center",fontsize=6,color="#555555")
labels=[]
for r,yy in zip(rows,ys):
    a1,b1=iv(r["breakpoint1"]); a2,b2=iv(r["breakpoint2"]); m1=(a1+b1)/2; m2=(a2+b2)/2
    ax.add_patch(Rectangle((1,yy-.3),GEN-1,.6,color=BL))
    ax.add_patch(Rectangle((m1,yy-.3),m2-m1,.6,color=OR))
    for a,b in ((a1,b1),(a2,b2)):
        ax.add_patch(Rectangle((a,yy-.3),max(b-a,6),.6,fill=False,hatch="////",edgecolor="black",linewidth=0))
    labels.append(f'{r["accession"]}  (in: {r["inside_parent"]} | out: {r["outside_parent"]})  fit {float(r["mosaic_fit_%"]):.0f}%')
ax.set_yticks(ys); ax.set_yticklabels(labels,fontsize=7)
for t,r in zip(ax.get_yticklabels(),rows):
    if float(r["mosaic_fit_%"])<80: t.set_color("#B00020")
for g,v in gy.items(): ax.text(GEN+30,sum(v)/len(v),SHORT[g],va="center",fontsize=8,fontweight="bold")
ax.set_xlim(-30,2480); ax.set_ylim(y+0.2,-1.1); ax.set_xticks(range(0,2001,250))
ax.set_xlabel("Genome position (nt, HM162346.1 numbering; circular genome, grey = intergenic regions)",loc="left")
ax.legend(handles=[Patch(color=OR,label="arc from 'inside' parent"),Patch(color=BL,label="remainder (wraps through origin) from other parent"),
    Patch(facecolor="white",edgecolor="black",hatch="////",label="breakpoint interval (between flanking informative sites)")],
    loc="lower center",bbox_to_anchor=(0.45,1.0),ncol=3,fontsize=8,frameon=False)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.text(0.01,0.005,"Parents are named per genome; colours show arc role, not lineage. Red labels: mosaic fit <80% of informative sites. "
    "Origin-side breakpoints at the left edge lie in an uninformative stretch through the origin (3SEQ places it in the last 30 alignment columns).",fontsize=6.5,color="#444444")
plt.tight_layout(rect=(0,0.02,1,1)); plt.savefig("fig_breakpoint_map_v3.png",dpi=300); plt.savefig("fig_breakpoint_map_v3.pdf")
print("written fig_breakpoint_map_v3.png / .pdf")
