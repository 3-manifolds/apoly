"""
As per Corollary 1.3 of http://arxiv.org/abs/math/0303017, the
Alexander polynomial for an L-space knot in S^3 is highly special.
"""

from sage.all import *

def one_possible_alex(n):
    assert min(n) > 0 and n == sorted(n) and len(set(n)) == len(n)
    t = PolynomialRing(QQ, 't').gen()
    k, m = len(n), max(n)
    terms = [(-1)**(k - j + 1) * (t**(m - n[j]) + t**(n[j] + m)) for j in range(len(n))]
    return (-1)**k * t**m+ sum(terms)

def possible_alex(N):
    for P in Partitions(N):
        for S in Permutations(P, len(P)):
            ans = [S[0]]
            for s in S[1:]:
                ans.append(ans[-1] + s)              
            yield one_possible_alex(ans)

def test(N):
    for f in possible_alex(N):
        bad_roots = [r for r, m in f.roots(CC) if m > 1]
        if bad_roots:
            print f
            for r in bad_roots:
                print abs(r), r
            good_roots = [r for r, m in f.roots(CC) if m == 1 and abs(abs(r) - 1) < 1e-10]
            print good_roots, '\n'
            


    



    
