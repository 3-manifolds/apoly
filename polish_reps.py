from sage.all import *
import sys, os, re, tempfile, random, string
import snappy, snappy.snap, euler
from shapes import polished_tetrahedra_shapes

def random_word(letters, N):
    return ''.join( [random.choice(letters) for i in range(N)] )

def inverse_word(word):
    return word.swapcase()[::-1]

def words_in_Fn(gens, n):
    next_letter = dict()
    sym_gens = gens + [g.swapcase() for g in gens]
    for g in sym_gens:
        next_letter[g] = [h for h in sym_gens if h != g.swapcase()]
    if n == 1:
        return sym_gens
    else:
        words = words_in_Fn(gens, n - 1)
        ans = []
        for word in words:
            ans += [word + g for g in next_letter[word[-1]] if len(word) == n - 1]
        return words + ans

def is_lex_first_in_conjugacy_class(word):
    if word[0] == word[-1].swapcase():
        return False
    for i in range(len(word)):
        other = word[i:] + word[:i]
        if other < word or other.swapcase() < word:
            return False
    return True
    
def conjugacy_classes_in_Fn(gens, n):
    return [word for word in words_in_Fn(gens, n) if is_lex_first_in_conjugacy_class(word)]

def SL2C_inverse(A):
    return matrix([[A[1,1], -A[0,1]], [-A[1,0], A[0, 0]]])

def GL2C_inverse(A):
    return (1/det(A))*matrix([[A[1,1], -A[0,1]], [-A[1,0], A[0, 0]]])

def apply_representation(word, gen_images):
    gens = string.ascii_lowercase[:len(gen_images)]
    rho = dict([(g, gen_images[i]) for i, g in enumerate(gens)] +
               [(g.upper(), SL2C_inverse(gen_images[i])) for i, g in enumerate(gens)])
    return prod( [rho[g] for g in word], Id2)

Id2 = MatrixSpace(ZZ, 2)(1)


class CheckRepresentationFailed(Exception):
    pass
    
class PSL2CRepOf3ManifoldGroup:
    """
    Throughout precision is in bits.
    """
    def __init__(self, manifold,
                 target_meridian_holonomy_arg = None,
                 rough_shapes=None,
                 precision=100,
                 fundamental_group_args=tuple() ):
        self.precision = precision
        self.manifold, self.rough_shapes = manifold.copy(), rough_shapes
        if target_meridian_holonomy_arg is None:
            RR = RealField()
            target_meridian_holonomy_arg = RR(
                manifold.cusp_info('holonomies')[0][0].imag) 
        self.target_meridian_holonomy_arg = target_meridian_holonomy_arg
        self.fundamental_group_args = fundamental_group_args
        self._cache = {}
        if rough_shapes != None:
            self.manifold.set_tetrahedra_shapes(rough_shapes, rough_shapes) 

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
                G = XXXpolished_holonomy(self.manifold,
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

def real_part_of_matricies_with_error(matrices):
    real_with_errors = [real_part_of_matrix_with_error(A) for A in matrices]
    return [r for r,e in real_with_errors], max(e for r, e in real_with_errors)

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

def right_kernel_two_by_two(A):
    """
    For a 2x2 matrix A over an approximate field like RR or CC find an
    element in the right kernel.
    """
    CC = A.base_ring()
    prec = CC.precision()
    epsilon = (RealField()(2.0))**(-0.8*prec)
    assert abs(A.determinant()) < epsilon
    a, b = A[0]
    v = vector(CC, [1, -a/b]) if abs(b) > abs(a) else vector(CC, [-b/a, 1])
    assert (A*v).norm() < epsilon
    return v/v.norm()
    
def eigenvectors(A):
    """
    Returns the two eigenvectors of a loxodromic matrix A
    """
    CC = A.base_ring()
    return [right_kernel_two_by_two(A-eigval) for eigval in A.charpoly().roots(CC, False)]
    
def eigenbasis(A, B):
    """
    Given loxodromic matrices A and B, return a basis of C^2 consisting of
    one eigenvector from each. 
    """
    basis = [ (a, b) for a in eigenvectors(A) for b in eigenvectors(B) ]
    return matrix(min(basis, key=lambda (a,b) : abs(a*b))).transpose()

def conjugator_into_PSL2R(A, B):
    """
    Given loxodromic matrices A and B which lie in a common conjugate of
    PSL(2, R), return a matrix C so that C^(-1)*A*C and C^(-1)*B*C are in
    PSL(2, R) itself.
    """
    C = eigenbasis(A, B)
    AA = GL2C_inverse(C)*A*C
    return C * matrix(A.base_ring(), [[1, 0], [0, 1/AA[0,1]]])

def conjugate_into_PSL2R(rho, max_error, depth=5):
    gens = rho.generators()
    new_mats, error = real_part_of_matricies_with_error(rho(g) for g in gens)
    if error < max_error:
        return new_mats

    for word in conjugacy_classes_in_Fn(gens, depth):
        U = rho(word)
        if abs(U.trace()) > 2.0001:
            conjugates = [ rho(g)*U*rho(g.upper()) for g in gens ]
            V = max(conjugates, key=lambda M: (U - M).norm())
            C =  conjugator_into_PSL2R(U, V)
            new_mats = [GL2C_inverse(C) * rho(g) * C for g in gens]
            final_mats, error = real_part_of_matricies_with_error(new_mats)
            assert error < max_error
            return final_mats

    raise ValueError("Couldn't conjugate into PSL(2, R)")

class PSL2RRepOf3ManifoldGroup(PSL2CRepOf3ManifoldGroup):
    def __init__(self, rep_or_manifold,
                 target_meridian_holonomy_arg=None,
                 rough_shapes=None,
                 precision=None,
                 fundamental_group_args=tuple()):
        if isinstance(rep_or_manifold, PSL2CRepOf3ManifoldGroup):
            rep = rep_or_manifold
        else:
           rep = PSL2CRepOf3ManifoldGroup(
               rep_or_manifold,
               target_meridian_holonomy_arg,
               rough_shapes,
               precision,
               fundamental_group_args)
        self.manifold = rep.manifold
        self.target_meridian_holonomy_arg = rep.target_meridian_holonomy_arg
        self.rough_shapes =  rep.rough_shapes
        self.precision = rep.precision
        self.fundamental_group_args = rep.fundamental_group_args
        self._cache = {}

    def polished_holonomy(self, precision=None):
        self._update_precision(precision)
        precision = self.precision
        if precision == None:
            raise ValueError, "Need to have a nontrivial precision set"
        mangled = "polished_holonomy_%s" % precision
        if not self._cache.has_key(mangled):
            epsilon = RR(2.0)**(-0.8*precision)
            G = XXXpolished_holonomy(self.manifold, self.target_meridian_holonomy_arg,
                                     precision,
                                     fundamental_group_args=self.fundamental_group_args,
                                     lift_to_SL2=False,
                                     ignore_solution_type=True)
            new_mats = conjugate_into_PSL2R(G, epsilon)
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
        shapes = "[" + ",".join(["%s" % z for z in self.rough_shapes]) + "]"
        traces = "[" + ",".join(["%s" % z for z in self.trace_field_generators()]) + "]"
        return "<%s" % self.manifold + ": " + traces + ">"


from snappy.snap import  generators
from snappy.snap.polished_reps import (initial_tet_ideal_vertices,
                                       reconstruct_representation,
                                       clean_matrix,
                                       ManifoldGroup)

def XXXpolished_holonomy(M, target_meridian_holonomy_arg,
                         bits_prec=100,
                         fundamental_group_args = [],
                         lift_to_SL2 = True,
                         ignore_solution_type=False,
                         dec_prec=None):
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

import euler

def lift_on_cusped_manifold(rho):
    rels = rho.relators()[:-1]
    euler_cocycle = [euler.euler_cocycle_of_relation(rho, R) for R in rels]
    D = rho.coboundary_1_matrix()[:-1]
    M = matrix(ZZ, [euler_cocycle] + D.columns())
    k = M.left_kernel().basis()[0]
    assert k[0] == 1
    shifts = (-k)[1:]
    good_lifts = [euler.PSL2RtildeElement(rho(g), s)
                  for g, s in zip(rho.generators(), shifts)]
    rho_til= euler.LiftedFreeGroupRep(rho, good_lifts)
    return rho_til

