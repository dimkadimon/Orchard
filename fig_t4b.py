import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
d=json.load(open('atleast4_N60_certificate.json'))
A,B,C,D=d['A'],d['B'],d['C'],d['D']
xs=list(range(13)); ys=list(range(13))
Cs=set(c-240 for c in C); Ds=set(e-120 for e in D)
ok=[(x,y) for x in xs for y in ys if (x+y) in Cs and (x-y) in Ds]
print("transversals:",len(ok)); assert len(ok)==129
fig,ax=plt.subplots(figsize=(5.0,5.0),dpi=200)
# strips
ax.add_patch(Rectangle((-0.5,-0.5),13,13,fc='#f7f7f9',ec='none',zorder=0))
for x in xs:
    for y in ys:
        good=(x+y) in Cs and (x-y) in Ds
        ax.add_patch(Rectangle((x-0.45,y-0.45),0.9,0.9,
            fc=('#1f77b4' if good else '#ffffff'),ec='#cccccc',lw=0.4,zorder=1))
for x,y in ok:
    ax.plot(x,y,'o',ms=3.0,color='white',zorder=2)
ax.set_xlim(-0.7,12.7); ax.set_ylim(-0.7,12.7)
ax.set_xticks(xs); ax.set_yticks(ys); ax.tick_params(labelsize=7,length=2)
ax.set_aspect('equal')
ax.set_xlabel(r'$x$  (index on axis $A$)',fontsize=9)
ax.set_ylabel(r'$y$  (index on axis $B$)',fontsize=9)
for s in ax.spines.values(): s.set_color('#999999'); s.set_linewidth(0.6)
ax.grid(False)
ax.set_title(r'$n=60$: the $129$ transversals are the pairs $(x,y)$'+'\n'+r'with $x+y\in C$ and $x-y\in D$',fontsize=9.5,pad=7)
fig.tight_layout(pad=0.5)
fig.savefig('paper/figures/fig_T4_60.pdf',bbox_inches='tight')
fig.savefig('paper/figures/fig_T4_60.png',bbox_inches='tight',dpi=200)
print("replaced T4 figure")
