from sage.all import *
from data.lspace import names, L_space_slopes, genus
import snappy


manifolds = [snappy.Manifold(name) for name in names]

def compute_genus(M):
    g0 = snappy.snap.alexander_polynomial(M).degree()/2
    g1 = (snappy.snap.hyperbolic_torsion(M, bits_prec=1000).degree()/2 + 1)/2
    assert g0 == g1
    return g0

def normalize(slope):
    return -slope if slope[0] < 0 else slope

class Cone:
    def __init__(self, v, w):
        self.a, self.b = v[0], v[1]
        self.c, self.d = w[0], w[1]
        if self.det() != 1:
            self.a, self.b = w[0], w[1]
            self.c, self.d = v[0], v[1]
        #if self.det() != 1:
        #    raise ValueError('Input vectors are not a basis for Z^2')
         
    def det(self):
        return self.a*self.d - self.b*self.c
    
    def _contains(self, x, y):
        u = self.d*x - self.c*y
        v = -self.b*x + self.a*y
        return u*v >= 0

    def gens(self):
        return vector( (self.a, self.b) ) , vector( (self.c, self.d) )
    
    def __contains__(self, s):
        return self._contains(s[0], s[1])

    def __repr__(self):
        return "<Cone (%d, %d) (%d, %d)>" % (self.a, self.b, self.c, self.d)


def compute_L_space_range(name):
    M = snappy.Manifold(name)
    K = snappy.CensusKnots.identify(M)
    A = M.is_isometric_to(K, True)[0].cusp_maps()[0]
    A = matrix(ZZ,  [[A[0,0], A[0,1]], [A[1,0], A[1,1]]])
    Ainv = A**(-1)
    try:
        g = genus[name]
        L_slopes = [vector(s) for s in L_space_slopes[name]]
        L_slopes_in_K = [normalize(A*s) for s in L_slopes]
        s = 1 if min(s[1] for s in L_slopes_in_K) >= 0 else -1
        knot_meridian =  vector(ZZ, (1, 0))
        l_space_edge = vector(ZZ, (2*g - 1, s))
    except KeyError:
        g = compute_genus(M)
        knot_meridian =  vector(ZZ, (2*g - 1, 1))
        l_space_edge = vector(ZZ, (2*g - 1, -1))
       
    C = Cone(Ainv*knot_meridian, Ainv*l_space_edge)
    #assert {slope in C for slope in L_space_slopes[name]} == {True}
    return C

