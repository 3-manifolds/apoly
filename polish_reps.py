from sage.all import *
import sys, os, re, tempfile, random, string
import snappy, snappy.snap, euler

def random_word(letters, N):
    return ''.join( [random.choice(letters) for i in range(N)] )

def SL2C_inverse(A):
    B = copy(A)
    B[0,0], B[1,1] = A[1,1], B[0,0]
    B[0,1], B[1,0] = -A[0,1], -B[1,0]
    return B

def apply_representation(word, gen_images):
    gens = string.ascii_lowercase[:len(gen_images)]
    rho = dict([(g, gen_images[i]) for i, g in enumerate(gens)] +
               [(g.upper(), SL2C_inverse(gen_images[i])) for i, g in enumerate(gens)])
    return prod( [rho[g] for g in word], Id2)

Id2 = MatrixSpace(ZZ, 2)(1)


class CheckRepresentationFailed(Exception):
    pass
    
class PSL2CRepOf3ManifoldGroup:
    def __init__(self, manifold, rough_shapes=None, precision=100):
        self.manifold, self.rough_shapes, self.precision = manifold.copy(), rough_shapes, precision
        self._cache = {}
        if rough_shapes != None:
            self.manifold.set_tetrahedra_shapes(rough_shapes, fillings = manifold.cusp_info('filling'))

    def __repr__(self):
        return "<%s" % self.manifold + ": [" + ",".join(["%.3f" % z for z in self.rough_shapes]) + "]>"

    def _update_precision(self, precision):
        if precision != None:
            self.precision = precision
            
        
    def polished_holonomy(self, precision=None):
        self._update_precision(precision)
        precision = self.precision
        mangled = "polished_holonomy_%s" % precision
        if not self._cache.has_key(mangled):
            if precision == None:
                G = self.manifold.fundamental_group()
            else:
                G = snappy.snap.polished_holonomy(self.manifold, bits_prec=precision,
                                                   lift_to_SL2=False, ignore_solution_type=True)
                if not G.check_representation() < RR(2.0)**(-0.8*precision):
                    raise CheckRepresentationFailed

            self._cache[mangled] = G
                    
        return self._cache[mangled]

    def trace_field_generators(self, precision=None):
        self._update_precision(precision)
        G = self.polished_holonomy()
        return G.trace_field_generators()

    def invariant_trace_field_generators(self, precision=None):
        self._update_precision(precision)
        G = self.polished_holonomy()
        return G.invariant_trace_field_generators()

    def has_real_traces(self, precision=None):
        self._update_precision(precision)
        real_precision = self.precision if self.precision else 15
        max_imaginary_part = max([ abs(tr.imag()) for tr in self.trace_field_generators()] )
        return  max_imaginary_part  <  RR(10.0)**(-0.5*real_precision)

    def appears_to_be_SU2_rep(self, trys=100, N = 20):
        G = self.polished_holonomy()
        for i in range(trys):
            w = random_word(G.generators(), N)
            d = abs(self(w).trace())
            if d > 2.1:
                return False

        return True

    def is_PSL2R_rep(self):
        rt = self.has_real_traces()
        not_su2 = not self.appears_to_be_SU2_rep()
        from_filling = self.really_comes_from_filling()
        return rt and not_su2 and from_filling

    def really_comes_from_filling(self):
        G = self.polished_holonomy()
        return G.check_representation() < RR(2.0)**(-0.8*self.precision)

    def generators(self):
        return self.polished_holonomy().generators()

    def relators(self):
        return self.polished_holonomy().relators()

    def name(self):
        return repr(self.manifold)

    def coboundary_1_matrix(self):
        gens, rels = self.generators(), self.relators()
        return matrix(ZZ, [[R.count(g) - R.count(g.swapcase()) for g in gens] for R in rels])

    def H2(self):
        """
        Computes H^2(G; Z) *assuming* d_3 : C_3 -> C_2 is the
        zero map. 
        """
        if not 'smith_form' in self._cache:
            self._cache['smith_form'] = self.coboundary_1_matrix().smith_form()

        D, U, V = self._cache['smith_form']
        ans = [d for d in D.diagonal() if d != 1]
        assert self.manifold.homology().coefficients == ans
        return ans

    def has_2_torsion_in_H2(self):
        H2 = self.H2()
        return len([c for c in H2 if c != 0 and c % 2 == 0]) > 0

    def class_in_H2(self, cocycle):
        self.H2()
        D, U, V = self._cache['smith_form']
        ans = []
        for c, d in zip(U*cocycle, D.diagonal()):
            if d != 1:
                a = c if d == 0 else c % d
                ans.append(a)
        return vector(ans)

    def matrix_field(self):
        return self.trace_field_generators()[0].parent()

    def __call__(self, word):
        return self.polished_holonomy().SL2C(word)


def real_part_of_matrix_with_error(A):
    RR = RealField(A.base_ring().precision())
    entries = A.list()
    real_parts, error = [x.real() for x in entries], max([abs(x.imag()) for x in entries])
    B = matrix(RR, real_parts, nrows=A.nrows(), ncols=A.ncols())
    if B.trace() < 0:
        B = -B
    return B, error

def normalize_vector(v):
    return v/v.norm()

def apply_matrix(mat, v):
    return normalize_vector(mat*v)

def swapped_dot(a, b):
    return -a[0]*b[1] + a[1]*b[0]

def orientation(a, b, c):
    return cmp( swapped_dot(a,b) * swapped_dot(b,c) * swapped_dot(c, a), 0)

def dist(a,b):
    return (a - b).norm()

class PSL2RRepOf3ManifoldGroup(PSL2CRepOf3ManifoldGroup):
    def __init__(self, rep_or_manifold, rough_shapes=None, precision=None):
        if isinstance(rep_or_manifold, PSL2CRepOf3ManifoldGroup):
            rep = rep_or_manifold
        else:
            rep = PSL2CRepOf3ManifoldGroup(rep_or_manifold, rough_shapes, precision)

        self.manifold, self.rough_shapes, self.precision = rep.manifold, rep.rough_shapes, rep.precision
        self._cache = {}

    def polished_holonomy(self, precision=None):
        self._update_precision(precision)
        precision = self.precision
        if precision == None:
            raise ValueError, "Need to have a nontrivial precision set"

        mangled = "polished_holonomy_%s" % precision
        if not self._cache.has_key(mangled):
            epsilon = RR(10.0)**(-0.8*precision)
            G = snappy.snap.polished_holonomy(self.manifold, precision,
                                               lift_to_SL2=False, ignore_solution_type=True)
            real_with_errors = [real_part_of_matrix_with_error(G.SL2C(g)) for g in G.generators()]
            if max( [e for r, e in real_with_errors] )  > epsilon:
                # This sometimes works:
                CC = G.SL2C('a').base_ring()
                i = CC.gen()
                A = matrix(CC, [[i,0],[0,1]])
                real_with_errors = [real_part_of_matrix_with_error(
                    A*G.SL2C(g)*(A**-1)) for g in G.generators()]
                assert [e for r,e in real_with_errors] < epsilon
                                                          
            new_mats = [r for r, e in real_with_errors]
            def rho(word):
                return apply_representation(word, new_mats)
            G.SL2C = rho
            if not G.check_representation() < epsilon:
                  raise CheckRepresentationFailed
            self._cache[mangled] = G
                    
        return self._cache[mangled]

    def thurston_class_of_relation(self, word, init_pt):
        """
        The Thurston Class is twice the Euler class.  Not sure WTF this means when there's
        2-torsion in H^2.  
        """
        n = len(word)
        b = normalize_vector(init_pt)
        points = [apply_matrix(self(word[:i]), b) for i in range(1, n)]
        error = min( [dist(b, p) for p in points] + [dist(points[i], points[i + 1]) for i in range(n - 2)])
        return sum( [ orientation(b, points[i], points[i+1]) for i in range(n - 2)]), error

    def thurston_class(self, init_pt = (2,-3)):
        init_pt = vector(self.matrix_field(), init_pt)
        ans = [self.thurston_class_of_relation(R, init_pt) for R in self.relators()]
        thurston, error = vector(ZZ, [x[0] for x in ans]), min([x[1] for x in ans])
        return self.class_in_H2(thurston), error

    def euler_class(self, double=False):
        rels = self.relators()
        e = vector(ZZ, [euler.euler_cocycle_of_relation(self, R) for R in rels])
        if double:
            e = 2*e
        return self.class_in_H2(e)

    def representation_lifts(self, precision=None):
        self._update_precision(precision)
        thurston, error = self.thurston_class()
        if thurston == 0:
            if not self.has_2_torsion_in_H2():
                return True
            else:
                return self.euler_class() == 0

        return False
                
    def __repr__(self):
        shapes = "[" + ",".join(["%.3f" % z for z in self.rough_shapes]) + "]"
        traces = "[" + ",".join(["%.3f" % z for z in self.trace_field_generators()]) + "]"
        return "<%s" % self.manifold + ": " + traces + ">"

    
