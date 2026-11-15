import numpy as np, sys
from fractions import Fraction
def vid(x,y,n1,n2,s):
    k=y//n2; x-=k*s; y-=k*n2; return (x%n1)+n1*y
for n1,n2,s in [(3,3,0),(3,3,1),(4,3,0)]:
    N=n1*n2; L=np.zeros((N,N))
    for y in range(n2):
        for x in range(n1):
            u=x+n1*y
            for v in (vid(x+1,y,n1,n2,s),vid(x,y+1,n1,n2,s)):
                L[u,u]+=1;L[v,v]+=1;L[u,v]-=1;L[v,u]-=1
    print(n1,n2,s,round(np.linalg.det(L[1:,1:])))
