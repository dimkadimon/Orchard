"""Publication figures for the six new algebraic 4-orchard records."""
import json, numpy as np, sympy as sp, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from algfield import parse_entry
from orchard import float_lines

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
    best=None
    for r,P in numeric_points(F,pts):
        if len(set(map(tuple,np.round(P,9))))<len(P): continue
        sizes=float_lines(P,tol=1e-7)
        if any(s>=5 for s in sizes): continue
        sc=len(sizes.get(4,[]))
        if best is None or sc>best[0]: best=(sc,r,P,sizes)
    assert best is not None, N
    sc,r,P,sizes=best
    print(f"N={N}: score={sc} (expected {expected}) root a={r:.10f} lines4={len(sizes.get(4,[]))} lines3={len(sizes.get(3,[]))}")
    assert sc==expected,(N,sc,expected)
    fig,ax=plt.subplots(figsize=(5.0,5.0),dpi=200)
    quads=sizes.get(4,[])
    cmap=plt.get_cmap('turbo')
    for k,q in enumerate(quads):
        Q=P[list(q)]
        ctr=Q.mean(axis=0); ang=np.arctan2(Q[:,1]-ctr[1],Q[:,0]-ctr[0])
        Q=Q[np.argsort(ang)]
        i0,i1=0,0; bd=-1
        for i in range(len(Q)):
            for j in range(i+1,len(Q)):
                d=np.linalg.norm(Q[i]-Q[j])
                if d>bd: bd=d;i0,i1=i,j
        A,B=Q[i0],Q[i1]
        ax.plot([A[0],B[0]],[A[1],B[1]],'-',color=cmap(k/max(1,len(quads)-1)),lw=0.55,alpha=0.75,zorder=1)
    ax.plot(P[:,0],P[:,1],'o',ms=3.0,mfc='#111111',mec='white',mew=0.4,zorder=3)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(f'$N={N}$,  $t_4^{{(5)}}\\geq {sc}$',fontsize=11,pad=6)
    ax.margins(0.06)
    fig.tight_layout(pad=0.4)
    fig.savefig(f'paper/figures/{fname}.pdf',bbox_inches='tight')
    fig.savefig(f'paper/figures/{fname}.png',bbox_inches='tight',dpi=200)
    plt.close(fig)
    return {'N':N,'score':sc,'root':r,'n4':len(quads),'n3':len(sizes.get(3,[]))}

if __name__=='__main__':
    info=[]
    for N,fn,exp in [(28,'fig_N28',42),(36,'fig_N36',63),(37,'fig_N37',63),
                     (40,'fig_N40',70),(48,'fig_N48',96),(49,'fig_N49',97)]:
        info.append(make(N,None,fn,exp))
    json.dump(info,open('paper/figures/figure_info.json','w'),indent=1)
    print("figures written")
