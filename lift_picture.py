from apoly import *
from polish_reps import *
from euler import *

def lift_on_cusped_manifold(rho):
#???    rels = rho.relators()[:-1]
    rels = rho.relators()
    euler_cocycle = [euler.euler_cocycle_of_relation(rho, R) for R in rels]
#???    D = rho.coboundary_1_matrix()[:-1]
    D = rho.coboundary_1_matrix()
    M = matrix(ZZ, [euler_cocycle] + D.columns())
    k = M.left_kernel().basis()[0]
    assert k[0] == 1
    shifts = (-k)[1:]
    good_lifts = [euler.PSL2RtildeElement(rho(g), s)
                  for g, s in zip(rho.generators(), shifts)]
    rho_til= euler.LiftedFreeGroupRep(rho, good_lifts)
    return rho_til

def elliptic_fixed_point(A):
    assert abs(A.trace()) < 2.0
    R = A.base_ring()
    C = R.complex_field()
    x = PolynomialRing(R, 'x').gen()
    a, b, c, d = A.list()
    p = c*x*x + (d - a)*x - b
    if p == 0:
        return CC.gen()
    return max(p.roots(CC, False), key=lambda z:z.imag())

def elliptic_rotation_angle(A):
    z = elliptic_fixed_point(A)
    
    a, b, c, d = A.list()
    derivative = 1/(c*z + d)**2
    pi = A.base_ring().pi()
    r = -derivative.argument()
    if r < 0:
        r = r + 2*pi
    return r/(2*pi)
    
def translation_amount(A_til):
    return elliptic_rotation_angle(A_til.A) + A_til.s

def rot(R, t, s):
    t = R.pi()*R(t)
    A = matrix(R, [[cos(t), -sin(t)], [sin(t), cos(t)]])
    return euler.PSL2RtildeElement(A, s)

def shift_of_central(A_til):
    assert A_til.is_central()
    return A_til.s

def in_SL2R(H, f, s):
    shape = H.T_fibers[f].shapes[s]
    ev = H.T_longitude_evs[s][f]
    if abs(1.0 - abs(ev)) > .00001:
        return False
    if H.in_SU2(shape):
        return False
    return True

def lifted_slope(M,  target_meridian_holonomy_arg, shapes):
    RR = RealField()
    target = RR(target_meridian_holonomy_arg)
    rho = PSL2CRepOf3ManifoldGroup(M, target, rough_shapes=shapes, precision=1000)
    assert rho.polished_holonomy().check_representation() < 1.0e-100
    rho_real = PSL2RRepOf3ManifoldGroup(rho)
    meridian, longitude = rho.polished_holonomy().peripheral_curves()[0]
    rho_tilde = lift_on_cusped_manifold(rho_real)
    return ( -translation_amount(rho_tilde(longitude)) /
             translation_amount(rho_tilde(meridian)) )

class SL2RLifter:
    def __init__(self, V):
        self.holonomizer = H = V.holonomizer
        self.degree = H.degree
        self.find_shapes()
        print 'lifting reps'
        self.find_reps()
        print 'computing translations'
        self.find_translation_arcs()

    def find_shapes(self):
        self.SL2R_arcs = []
        H = self.holonomizer
        current_arc = None
        for s in range(self.degree):
            saving = False
            for n in range(128):
                try:
                    point_is_good = in_SL2R(H, n, s)
                except:
                    point_is_good = False
                if not saving and point_is_good:
                    current_arc = [ ((s,n), H.T_fibers[n].shapes[s]) ]
                    saving = True
                elif saving and point_is_good:
                    current_arc.append( ((s,n), H.T_fibers[n].shapes[s]) )
                elif saving and not point_is_good:
                    if len(current_arc) > 1:
                        self.SL2R_arcs.append(current_arc)
                    #elif len(current_arc) == 1:
                    #    print 'skipping', n-1, s
                    saving = False

    def find_reps(self):
        self.SL2R_rep_arcs = []
        RR = RealField(1000)
        for arc in self.SL2R_arcs:
            reps = []
            for sn,  S in arc:
                s, n = sn
                target1 = log(self.holonomizer.T_fibers[n].H_meridian).imag()
                target = -2*RR(pi)*RR(n)/RR(128)
                rho = PSL2RRepOf3ManifoldGroup(
                    self.holonomizer.manifold,
                    target,
                    S,
                    precision=1000,
                    fundamental_group_args = [True, False, True])
                if rho.polished_holonomy().check_representation() < 1.0e-100:
                    reps.append( (sn, rho) )
            if len(reps) > 1: 
                self.SL2R_rep_arcs.append(reps)

    def find_translation_arcs(self):
        self.translation_arcs = []
        self.translation_dict = {}
        for arc in self.SL2R_rep_arcs:
            #print len(arc)
            translations = []
            for sn, rho in arc:
                meridian, longitude = rho.polished_holonomy().peripheral_curves()[0]
                rho_til = lift_on_cusped_manifold(rho)
                try:
                    P = ( float(translation_amount(rho_til(meridian))),
                          float(translation_amount(rho_til(longitude))) )
                    translations.append(P)
                    self.translation_dict[sn] = P
                except AssertionError:
                    pass
            self.translation_arcs.append(translations)

    def show(self):
        plotlist = [ [complex(x,y) for x, y in arc] for arc in self.translation_arcs ]
        plot = Plot(plotlist)

    def show_slopes(self):
        M = self.holonomizer.manifold.copy()
        plotlist = []
        for arc in self.SL2R_arcs:
            slopes = []
            for sn,  S in arc:
                M.set_tetrahedra_shapes(S, S, [(0,0)])
                Hm, Hl = M.cusp_info('holonomies')[0]
                if Hm.imag != 0:
                    slopes.append(float(-Hl.imag/Hm.imag))
                elif len(slopes) > 1:
                    slopes.append(None)
            plotlist.append(slopes)
        plot = Plot(plotlist)


