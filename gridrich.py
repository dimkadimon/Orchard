"""Exact count of lines with >=k points in lattice point sets of various shapes.
Method: for every primitive direction v=(a,b) in a half-plane, the lines with
direction v are the level sets of c = b*x - a*y.  Count level sets with >=k points.
This is exact and runs in O(#directions * N)."""
import numpy as np, math, json, sys

def primitive_dirs(rmax):
    out=[]
    R=int(math.ceil(rmax))
    for a in range(0,R+1):
        for b in range(-R,R+1):
            if a==0 and b<=0: continue
            if a*a+b*b>rmax*rmax: continue
            if math.gcd(a,abs(b))!=1: continue
            out.append((a,b))
    return out

def count_rich(pts,k):
    x=pts[:,0].astype(np.int64); y=pts[:,1].astype(np.int64)
    rmax=max(np.abs(x).max(),np.abs(y).max())*2.0/(k-1)+1
    total=0
    for a,b in primitive_dirs(rmax):
        c=b*x-a*y
        c-=c.min()
        cnt=np.bincount(c)
        total+=int((cnt>=k).sum())
    return total

def disk(r):
    g=np.arange(-r,r+1)
    X,Y=np.meshgrid(g,g,indexing='ij')
    m=X*X+Y*Y<=r*r
    return np.stack([X[m],Y[m]],axis=1)

def square(s):
    g=np.arange(s)
    X,Y=np.meshgrid(g,g,indexing='ij')
    return np.stack([X.ravel(),Y.ravel()],axis=1)

def hexagon(r):
    # points of Z^2 inside the hexagon |x|<=r, |y|<=r, |x+y|<=r  (affine-regular hexagon)
    g=np.arange(-r,r+1)
    X,Y=np.meshgrid(g,g,indexing='ij')
    m=(np.abs(X)<=r)&(np.abs(Y)<=r)&(np.abs(X+Y)<=r)
    return np.stack([X[m],Y[m]],axis=1)

def triangle(s):
    g=np.arange(s+1)
    X,Y=np.meshgrid(g,g,indexing='ij')
    m=(X+Y<=s)
    return np.stack([X[m],Y[m]],axis=1)

if __name__=='__main__':
    k=4
    print('shape param   N      rich>=4    ratio*N^2 ->  c = count/N^2   (target 4/(9pi^2)=%.6f)'%(4/(9*math.pi**2)))
    res={}
    for name,fn,params in [('disk',disk,[4,6,8,12,16,24,32,48,64,96]),
                           ('square',square,[8,12,16,24,32,48,64,96,128]),
                           ('hexagon',hexagon,[6,9,12,18,24,36,48,72]),
                           ('triangle',triangle,[10,16,24,36,48,72,96])]:
        rows=[]
        for p in params:
            pts=fn(p); N=len(pts)
            c=count_rich(pts,k)
            rows.append((p,N,c,c/N**2))
            print('%-9s %4d %7d %9d   %.6f'%(name,p,N,c,c/N**2))
        res[name]=rows
    json.dump(res,open('gridrich_results.json','w'),indent=1)
