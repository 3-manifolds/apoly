from sage.all import vector, matrix, MatrixSpace, ZZ, RR, RealField, prod, PolynomialRing
import sys, os, re, tempfile, random, string
from shapes import polished_tetrahedra_shapes

def SL2C_inverse(A):
    return matrix([[A[1,1], -A[0,1]], [-A[1,0], A[0, 0]]])

def GL2C_inverse(A):
    return (1/A.det())*matrix([[A[1,1], -A[0,1]], [-A[1,0], A[0, 0]]])

def apply_representation(word, gen_images):
    gens = string.ascii_lowercase[:len(gen_images)]
    rho = dict([(g, gen_images[i]) for i, g in enumerate(gens)] +
               [(g.upper(), SL2C_inverse(gen_images[i])) for i, g in enumerate(gens)])
    return prod( [rho[g] for g in word], Id2)

Id2 = MatrixSpace(ZZ, 2)(1)

def polished_holonomy(M, target_meridian_holonomy_arg,
                         bits_prec=100,
                         fundamental_group_args = [],
                         lift_to_SL2 = True,
                         ignore_solution_type=False,
                         dec_prec=None):

    from snappy.snap import  generators
    from snappy.snap.polished_reps import (initial_tet_ideal_vertices,
                                       reconstruct_representation,
                                       clean_matrix,
                                       ManifoldGroup)
    if dec_prec:
        bits_prec = None
        error = ZZ(10)**(-dec_prec*0.8)
    else:
        error = ZZ(2)**(-bits_prec*0.8)
    shapes = polished_tetrahedra_shapes(M, target_meridian_holonomy_arg, bits_prec=bits_prec, dec_prec=dec_prec)
    G = M.fundamental_group(*fundamental_group_args)
    N = generators.SnapPy_to_Mcomplex(M, shapes)
    init_tet_vertices = initial_tet_ideal_vertices(N)
    generators.visit_tetrahedra(N, init_tet_vertices)
    mats = generators.compute_matrices(N)
    gen_mats = [clean_matrix(A, error=error) for A in reconstruct_representation(G, mats)]
    PG = ManifoldGroup(G.generators(), G.relators(), G.peripheral_curves(), gen_mats)
    if lift_to_SL2:
        PG.lift_to_SL2C()
    else:
        assert PG.is_projective_representation()

    return PG


class CheckRepresentationFailed(Exception):
    pass
    
class PSL2CRepOf3ManifoldGroup:
    """
    Throughout precision is in bits.
    """
    def __init__(self, manifold,
                 target_meridian_holonomy_arg=None,
                 rough_shapes=None,
                 precision=100,
                 fundamental_group_args=tuple() ):
        self.precision = precision
        self.manifold, self.rough_shapes = manifold.copy(), rough_shapes
        if rough_shapes != None:
            self.manifold.set_tetrahedra_shapes(rough_shapes, rough_shapes) 
        if target_meridian_holonomy_arg is None:
            CC = ComplexField()
            holonomy = CC(complex(manifold.cusp_info('holonomies')[0][0]))
            target_meridian_holonomy_arg = holonomy.imag()
        self.target_meridian_holonomy_arg = target_meridian_holonomy_arg
        self.fundamental_group_args = fundamental_group_args
        self._cache = {}

    def __repr__(self):
        return "<%s" % self.manifold + ": [" + ",".join(["%s" % z for z in self.rough_shapes]) + "]>"

    def _update_precision(self, precision):
        if precision != None:
            self.precision = precision
        
    def polished_holonomy(self, precision=None):
        self._update_precision(precision)
        precision = self.precision
        mangled = "polished_holonomy_%s" % precision
        if not self._cache.has_key(mangled):
            if precision == None:
                G = self.manifold.fundamental_group(*self.fundamental_group_args)
            else:
                G = polished_holonomy(self.manifold,
                                         self.target_meridian_holonomy_arg,
                                         bits_prec=precision,
                                         fundamental_group_args=self.fundamental_group_args,
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
        return  max_imaginary_part  <  RR(2.0)**(-0.5*real_precision)

    def appears_to_be_SU2_rep(self, depth=5, trys=50, rand_length = 20):
        G = self.polished_holonomy()
        gens = G.generators()
        words = conjugacy_classes_in_Fn(gens, depth)
        words += [random_word(gens, rand_length) for i in range(trys)]
        for w in words:
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

    def peripheral_curves(self):
        return self.manifold.fundamental_group().peripheral_curves()

    def meridian(self):
        return self.peripheral_curves()[0][0]
    
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
