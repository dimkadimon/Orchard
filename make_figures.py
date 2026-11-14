"""Publication figures for the six new algebraic 4-orchard records."""
import json, numpy as np, sympy as sp, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from algfield import parse_entry
from orchard import analyze
from scipy.spatial.distance import pdist

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.6})

def numeric_points(F,pts):
    x=sp.Symbol('x')
    poly=sum(int(c)*x**i for i,c in enumerate(F.p))
    roots=[complex(r) for r in sp.Poly(poly,x).all_roots() if abs(sp.im(r))<1e-30]
    roots=sorted(set(float(sp.re(r)) for r in roots))
    out=[]
    for r in roots:
        P=[]
        for (u,v) in pts:
            fu=sum(float(c)*r**i for i,c in enumerate(u))
            fv=sum(float(c)*r**i for i,c in enumerate(v))
            P.append((fu,fv))
        out.append((r,np.array(P)))
    return out

def make(N,label,fname,expected):
    text=open(f'entries/N{N}.txt').read()
    F,pts=parse_entry(text)
    # Determine incidence lists exactly, then pick a real embedding with the
    # greatest minimum point separation relative to the diagram's span.
    exact = analyze(F, pts)
    assert exact['valid'] and exact['score'] == expected
    sizes = exact['sizes']
    options=[]
    for r,P in numeric_points(F,pts):
        separation = float(np.min(pdist(P)) / np.max(np.ptp(P, axis=0)))
        options.append((separation,r,P))
    separation,r,P=max(options,key=lambda v:v[0])
    sc=exact['score']
    print(f"N={N}: exact score={sc}, a={r:.10f}, relative minimum separation={separation:.4f}")
    fig,ax=plt.subplots(figsize=(5.2,5.2),dpi=250)
    quads=sizes.get(4,[])
    for q in quads:
        Q=P[list(q)]
        # The farthest pair is the visible segment spanned by its four points.
        i,j=max(((i,j) for i in range(4) for j in range(i+1,4)),
                key=lambda ij:np.linalg.norm(Q[ij[0]]-Q[ij[1]]))
        A,B=Q[i],Q[j]
        ax.plot([A[0],B[0]],[A[1],B[1]],'-',color='#69808c',
                lw=0.50,alpha=0.55,zorder=1)
    # A white halo prevents a dark or crowded line from swallowing a point.
    ax.scatter(P[:,0],P[:,1],s=56,facecolors='white',edgecolors='white',linewidths=1.2,zorder=3)
    ax.scatter(P[:,0],P[:,1],s=27,facecolors='#101820',edgecolors='none',zorder=4)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(f'$N={N}$,  $s_4\\geq {sc}$',fontsize=13,pad=8)
    ax.margins(0.06)
    fig.tight_layout(pad=0.4)
    fig.savefig(f'paper/figures/{fname}.pdf',bbox_inches='tight')
    fig.savefig(f'paper/figures/{fname}.png',bbox_inches='tight',dpi=200)
    plt.close(fig)
    return {'N':N,'score':sc,'root':r,'n4':len(quads),'n3':len(sizes.get(3,[]))}

if __name__=='__main__':
    info=[]
    # N=37 and N=49 are deliberately not drawn: each is the preceding set plus one
    # point, so the diagram adds nothing. They are also out of the paper.
    for N,fn,exp in [(28,'fig_N28',42),(36,'fig_N36',63),
                     (40,'fig_N40',70),(48,'fig_N48',96)]:
        info.append(make(N,None,fn,exp))
    json.dump(info,open('paper/figures/figure_info.json','w'),indent=1)
    print("figures written")
