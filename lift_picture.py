from apoly import Plot
from sage.all import RealField, ComplexField, log, pi
from manifold_reps.real_reps import PSL2RRepOf3ManifoldGroup, translation_amount

def in_SL2R(H, f, s):
    shape = H.T_fibers[f].shapes[s]
    ev = H.T_longitude_evs[s][f][1]
    if abs(1.0 - abs(ev)) > .00001:
        return False
    if H.in_SU2(shape):
        return False
    return True

class SL2RLifter:
    def __init__(self, V):
        self.holonomizer = H = V.holonomizer
        self.degree = H.degree
        self.manifold_name = V.manifold_name
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
                rho_til = rho.lift_on_cusped_manifold()
                if rho_til is None:
                    continue 
                try:
                    P = ( float(translation_amount(rho_til(meridian))),
                          float(translation_amount(rho_til(longitude))) )
                    if P[0] < 0:
                        P = (-P[0], -P[1])
                    translations.append(P)
                    self.translation_dict[sn] = P
                except AssertionError:
                    pass
            self.translation_arcs.append(translations)

    def show(self):
        plotlist = [ [complex(x,y) for x, y in arc] for arc in self.translation_arcs ]
        self.plot = Plot(plotlist, title=self.manifold_name)

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
        self.slope_plot = Plot(plotlist)

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

def check_slope(H, n, s):
    F = H.T_fibers[n]
    S = F.shapes[s]
    M = H.manifold.copy()
    target = log(F.H_meridian).imag()
    return float(lifted_slope(M, target, S))

def bisection(H, low, high, s, target_slope, epsilon=1.0e-8):
    CC = ComplexField()
    low_fiber = H.T_fibers[low]
    high_fiber = H.T_fibers[high]
    M = H.manifold
    F = H.fibrator
    assert check_slope(H,low,s) < target_slope < check_slope(H,high,s)
    print 'finding:', target_slope
    count = 0
    while count < 100:
        z = (low_fiber.H_meridian + high_fiber.H_meridian)/2
        target_holonomy = z/abs(z)
        target_holonomy_arg = CC(target_holonomy).log().imag()
        new_fiber = F.transport2(low_fiber, complex(target_holonomy))
        shapes = new_fiber.shapes[s]
        new_slope = lifted_slope(M, target_holonomy_arg, shapes)
        if abs(new_slope - target_slope) < epsilon:
            return new_fiber.shapes[s]
        if new_slope < target_slope:
            low_fiber = new_fiber
            print new_slope, 'too low'
        else:
            high_fiber = new_fiber
            print new_slope, 'too high'
        count += 1
        print count
    print 'limit exceeded'
    return new_fiber.shapes[s]
