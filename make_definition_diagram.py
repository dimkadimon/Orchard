#!/usr/bin/env python3
"""Reproduce the illustrative s4, t4, T4 line-census diagram."""
from pathlib import Path
from itertools import combinations
from collections import Counter
from math import gcd
from functools import reduce
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Both point sets are unions of the same two axes. Any non-axis line meets
# each axis at most once, so the displayed census is complete.
P7=[(x,0) for x in range(4)]+[(0,y) for y in range(1,4)]
P8=[(x,0) for x in range(5)]+[(0,y) for y in range(1,4)]
def histogram(P):
    lines={}
    for (x,y),(u,v) in combinations(P,2):
        a,b,c=y-v,u-x,x*v-u*y
        d=reduce(gcd,(a,b,c)); a,b,c=(a//d,b//d,c//d)
        if next(z for z in (a,b,c) if z)<0:a,b,c=-a,-b,-c
        lines.setdefault((a,b,c),set()).update(((x,y),(u,v)))
    return Counter(map(len,lines.values()))
assert histogram(P7)[4]==2 and max(histogram(P7))==4
assert histogram(P8)[4]==1 and histogram(P8)[5]==1 and max(histogram(P8))==5

blue='#21618c'; orange='#b05216'; ink='#102333'
fig, rows=plt.subplots(3,2,figsize=(7.2,7.3),gridspec_kw={'width_ratios':[1,1.9]})
for i,(ax,text) in enumerate(rows):
    P=P7 if i==0 else P8
    ax.set_aspect('equal');ax.set_xlim(-.45,4.55);ax.set_ylim(-.55,3.55);ax.axis('off')
    ax.plot([0,0],[-.15,3.15],color=blue,lw=2.6,zorder=1)
    hcol=orange if i==1 else blue
    ax.plot([-.15,3.15 if i==0 else 4.15],[0,0],color=hcol,lw=2.6,ls='--' if i==1 else '-',zorder=1)
    ax.scatter([p[0] for p in P],[p[1] for p in P],s=62,facecolor=ink,edgecolor='white',linewidth=.9,zorder=3)
    ax.text(-.17,3.23,'4',fontsize=10,color=blue,weight='bold')
    ax.text(2.9 if i==0 else 4.04,.20,'4' if i==0 else '5',fontsize=10,color=hcol,weight='bold')
    text.axis('off')
    if i==0:
        title=r'$s_4$: strict no-five variant'
        desc='7 points; both lines have exactly 4.\n2 counted lines; no line has 5 points.'
    elif i==1:
        title=r'$t_4$: exactly-four variant'
        desc='8 points; the 5-point line is allowed\nbut not counted. 1 counted line.'
    else:
        title=r'$T_4$: at-least-four variant'
        desc='Same 8 points; both the 4-point and\n5-point lines count. 2 counted lines.'
    text.text(.02,.70,title,fontsize=15,weight='bold',color=ink,va='center')
    text.text(.02,.40,desc,fontsize=11.5,color=ink,va='center',linespacing=1.35)
fig.text(.5,.015,'Solid blue: counted line     Dashed orange: not counted     Black dot: point',ha='center',fontsize=10.5,color=ink)
fig.subplots_adjust(left=.035,right=.995,top=.985,bottom=.085,hspace=.22,wspace=.04)
out=Path('paper/figures');out.mkdir(parents=True,exist_ok=True)
fig.savefig(out/'fig_definitions.pdf',bbox_inches='tight')
fig.savefig(out/'fig_definitions.png',dpi=220,bbox_inches='tight')
print('P7 complete histogram:',dict(sorted(histogram(P7).items())))
print('P8 complete histogram:',dict(sorted(histogram(P8).items())))
