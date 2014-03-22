"""
Using PHC to find solutions to the gluing equations for closed manifolds.
"""

import os, sys, re
import snappy, phc
from polish_reps import PSL2CRepOf3ManifoldGroup, CheckRepresentationFailed
from real_reps import PSL2RRepOf3ManifoldGroup

def clean_complex(z, epsilon=1e-20):
    r, i = abs(z.real), abs(z.imag)
    if r < epsilon and i < epsilon:
        ans = 0.0
    elif r < epsilon:
        ans = z.imag*1j
    elif i < epsilon:
        ans = z.real
    else:
        ans = z
    assert abs(z - ans) < epsilon
    return ans

class PHCGluingSolutionsOfClosed:
    def __init__(self, manifold):
        if isinstance(manifold, str):
            manifold = snappy.Manifold(manifold)
        else:
            manifold = manifold.copy()
        if True in manifold.cusp_info('complete?') or not manifold.is_orientable():
            raise ValueError("Manifold must be closed and orientable")

        manifold.set_peripheral_curves('fillings')
        self.manifold = manifold

        self.N = N = manifold.num_tetrahedra()
        variables = ['X%s'%n for n in range(N)] + ['Y%s'%n for n in range(N)]
        self.ring = phc.PolyRing(variables)
        self.equations = [self.rect_to_PHC(eqn) for eqn in
                          snappy.snap.shapes.enough_gluing_equations(manifold)]
        self.equations += ['X%s + Y%s - 1'%(n,n) for n in range(N)]
        self.system = phc.PHCSystem(self.ring,
                                    [phc.PHCPoly(self.ring, eqn) for eqn in self.equations])
        
    def raw_solutions(self, max_err=1e-6):
        ans = []
        for sol in self.system.solution_list():
            if sol.err < max_err:
                ans.append([clean_complex(z) for z in sol.point[:self.N]])
        return ans

    def solutions(self, working_prec=230):
        psl2Rtilde, psl2R, su2, rest = [], [], [], []
        for sol in self.raw_solutions():
            rho = PSL2CRepOf3ManifoldGroup(self.manifold,
                            target_meridian_holonomy_arg=0,
                            rough_shapes=sol)
            try:
                rho.polished_holonomy(working_prec)
            except:
                continue
            if rho.appears_to_be_SU2_rep():
                su2.append(sol)
            elif rho.is_PSL2R_rep():
                rho = PSL2RRepOf3ManifoldGroup(rho)
                if rho.representation_lifts():
                    psl2Rtilde.append(sol)
                else:
                    psl2R.append(sol)
            else:
                rest.append(sol)

        return psl2Rtilde, psl2R, su2, rest
                
    def rect_to_PHC(self, eqn):
        A, B, c = eqn
        left, right = ['1'], ['1']
        for n, a in enumerate(A):
            if a > 0:
                left += ['X%s'%n]*int(a)
            else:
                right += ['X%s'%n]*int(-a)
        for n, b in enumerate(B):
            if b > 0:
                left += ['Y%s'%n]*int(b)
            else:
                right += ['Y%s'%n]*int(-b)
        op = ' - ' if c == 1 else ' + '
        return '*'.join(left) + op + '*'.join(right)

if __name__=='__main__':
    M = snappy.Manifold('m004(1,2)')
    ans = PHCGluingSolutionsOfClosed(M).solutions()
    
