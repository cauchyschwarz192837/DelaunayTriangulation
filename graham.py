from primitives import *

def graham(P):
    P = sorted(P)

    U = P[:2]
    L = P[:2]
    for p in P[2:]:
        
        U.append(p)

        while len(U) >= 3 and ccw(U[-3], U[-2], U[-1]):
            del U[-2]
        
        L.append(p)

        while len(L) >= 3 and cw(L[-3], L[-2], L[-1]):
            del L[-2]

    L.extend(reversed(U[1:-1]))
    return L