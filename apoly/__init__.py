# -*- coding: utf-8 -*-

from random import random, randint
from subprocess import Popen, PIPE
import time, sys, os, Tkinter
import numpy
from numpy import array, matrix, ndarray
from numpy import dot, prod, diag, transpose, zeros, ones, eye
from numpy import log, exp, pi, sqrt, ceil
from numpy import dtype, take, arange, sum, where, vstack, argmin
from numpy import float64
from numpy.linalg import svd, norm, eig, solve, lstsq
from numpy.fft import ifft
try:
    from numpy.linalg import matrix_rank
except ImportError:
    def matrixrank(A,tol=1e-8):
        s = svd(A,compute_uv=0)
        return sum( where( s>tol, 1, 0 ) )
from phc import *
import snappy
snappy.SnapPy.matrix = matrix
snappy.SnapPyHP.matrix = matrix
from snappy import *
from spherogram.graphs import Graph
from point import PEPoint
#from plot import GnuplotPlot as Plot
from plot import MatplotPlot as Plot
global_linewidth=1

real_array = numpy.vectorize(float)

# Constants for Newton method
RESIDUAL_BOUND = 1.0E-14
STEPSIZE_BOUND = 1.0E-14

got_sage = False
try:
    from sage.all import polygen, polygens, QQ, CC
    M,L = polygens(QQ, 'M,L')
    t = polygen(QQ, 't')
    z = polygen(CC, 'z')
    got_sage = True
except:
    pass

def sage_2poly(dict, ring=None):
    if not got_sage:
        raise RuntimeError, 'SagePolynomials only exist when running in Sage.'
    if ring is None:
        ring = QQ['x','y']
    x, y = ring.gens()[:2]
    result = ring(0)
    for m, n in dict.keys():
        result += dict[(m,n)]*x**m*y**n
    return result
    
DTYPE = dtype('c16')

class GluingSystem:
    """
    The system of gluing equations for a manifold, with specified
    meridian holonomy.  If the manifold has n tetrahedra, we use the
    first n-1 edge equations, together with the equation for meridian
    holonomy.  The left hand side of each equation is a Laurent monomial
    in z_i and (1-z_i), where z_i are the shape paremeters.  The
    right hand side of the system is [1,...,1,Hm] where Hm is the
    the meridian holonomy (i.e. the first eigenvalue squared).
    """
    def __init__(self, manifold):
        self.manifold = manifold.copy()
        eqns = manifold.gluing_equations('rect')
        # drop the last edge equation
        self.glunomials = [Glunomial(A, B, c) for A, B, c in eqns[:-3]]
        self.M_nomial, self.L_nomial = [Glunomial(A, B, c) for A,B,c in eqns[-2:]]
        self.glunomials.append(self.M_nomial)

    def __repr__(self):
        return '\n'.join([str(G) for G in self.glunomials])

    def __call__(self, Z):
        return array([G(Z) for G in self.glunomials])

    def __len__(self):
        return len(self.glunomials)
    
    def jacobian(self, Z):
        return matrix([G.gradient(Z) for G in self.glunomials])

    def M_holonomy(self, Z):
        return complex(self.M_nomial(Z))

    def L_holonomy(self, Z):
        return complex(self.L_nomial(Z))
    
    def slope(self, Z, digits=10):
        """
        Use LLL to guess which peripheral element is being killed by
        the rep corresponding to shapes Z.  We look for small integers
        p, q and n such that p*phi + q*theta = 2*pi*n where phi and
        theta are the arguments of the meridian and longitude holonomies.
        """
        Hm, Hl = self.M_holonomy(Z), self.L_holonomy(Z)
        phi, theta = log(Hm).imag, log(Hl).imag
        N = 10**digits
        A = pari.matrix(4, 3,
                        [1, 0, 0,
                         0, 1, 0,
                         0, 0, 1,
                         phi*N, theta*N, -2*pi*N])
        L = A.qflll()
        return L[0][:2]

    def condition(self, Z):
        """
        Return the condition numbers of the Jacobians for the defining
        equations of the gluing curve and of the entire system
        at the point Z.
        """
        U, D, V = svd(self.jacobian(Z)[:-1])
        curve = D[0]/D[-1]
        U, D, V = svd(self.jacobian(Z))
        system = D[0]/D[-1]
        return curve, system
    
    def newton_step(self, Z, M_target):
        """
        Do one iteration of Newton's method, starting at Z and aiming
        to solve G(z) = (1,1,...,M_target).  Returns a triple:
        Z', step_size, residual.  Solves the linear system by 
        LU factorization (not great for nearly singular systems).
        """
        J = self.jacobian(Z)
        target = ones(len(self), dtype=DTYPE)
        target[-1] = M_target
        dZ = solve(J, target - self(Z))
        step_size = norm(dZ)
        Zn = Z + dZ
        return Zn, step_size, max(abs(target - self(Zn)))
    
    def newton_step_ls(self, Z, M_target):
        """
        Do one iteration of Newton's method, starting at Z and aiming
        to solve G(z) = (1,1,...,M_target). Returns a pair:
        dZ, (1,1,...,M_target).  Finds a least squares approcimate
        solution to the linear system.  This method is stable with nearly
        singular systems.
        """
        J = self.jacobian(Z)
        target = ones(len(self), dtype=DTYPE)
        target[-1] = M_target
        error = target - self(Z)
        dZ, residues, rank, sing = lstsq(J, error)
        return dZ, target

    def newton1(self, Z, M_target, debug=False):
        """
        Simple version of Newton's method.  Uses the LU decomposition to
        solve the linear system.  Does not adjust step sizes.
        The iteration is terminated if:
          * the residual does not decrease; or
          * the step size is smaller than 1.0E-15
          * more than 10 iterations have been attempted
        """
        prev_residual = step_size = 1.0E5
        prev_Z, count = Z, 1
        while True:
            Zn, step_size, residual = self.newton_step(prev_Z, M_target)
            if debug: print count, residual, step_size
            if residual > prev_residual:
                return prev_Z, prev_residual
            if (step_size < STEPSIZE_BOUND or
                residual < RESIDUAL_BOUND or
                count > 10):
                return Zn, residual
            prev_Z, prev_residual = Zn, residual
            count += 1

    def newton2(self, Z, M_target, debug=False):
        """
        Fancier version of Newton's method which uses a least squares
        solution to the linear system.  To avoid overshooting, step
        sizes are adjusted by an Armijo rule that successively halves
        the step size.
        The iteration is terminated if:
          * the residual does not decrease; or
          * the step size is smaller than 1.0E-15
          * more than 10 iterations have been attempted
        """
        prev_residual = step_size = 1.0E5
        prev_Z, count = Z, 1
        while True:
            dZ, target = self.newton_step_ls(prev_Z, M_target)
            t = 1.0
            for k in range(5):
                Zn = prev_Z + t*dZ
                residual = max(abs(target - self(Zn)))
                if residual < prev_residual:
                    break
                else:
                    t *= 0.5
            if debug: print 'scaled dZ by %s; residual: %s'%(t, residual)
            if residual > prev_residual:
                if debug: print 'Armijo failed with t=%s'%t
                return prev_Z, prev_residual
            step_size = norm(t*dZ)
            if (step_size < STEPSIZE_BOUND or
                residual < RESIDUAL_BOUND or
                count > 10):
                return Zn, residual
            prev_Z, prev_residual = Zn, residual
            count += 1

    def track(self, Z, M_target, dT=1.0, debug=False):
        """
        Track solutions of the gluing system starting at Z and
        ending at a solution where the meridian holonomy takes the
        specified value.  The path is subdivided into subpaths of
        length determined by dT, which may be further divided if
        convergence problems arise.
        """
        M_start = self(Z)[-1]
        delta = (M_target - M_start)
        T = 0.0
        Zn = Z
        if debug: print 'Z = %s; condition=%s'%(Z, [self.condition(x) for x in Z]) 
        # First we try the cheap and easy method
        while T < 1.0:
            T = T+dT
            target = M_start + T*delta
            Zn, residual = self.newton1(Zn, target)
        if residual < 1.0E-14: # What is a good threshold here?
            return Zn
        # If that fails, try taking baby steps.
        if debug: print 'Taking baby steps ...'
        success = 0
        T, dT = 0.0, 0.5*dT
        prev_Z = Z
        while T < 1.0:
            Tn = min(T+dT, 1.0)
            if debug: print 'trying T = %.17f'%Tn
            baby_target = M_start + Tn*delta
            Zn, residual = self.newton2(prev_Z, baby_target, debug=debug)
            if residual < 1.0E-12:
                prev_Z = Zn
                if success > 3:
                    success = 0
                    dT *= 2
                    if debug: print 'Track step increased to %.17f'%dT
                else:
                    success += 1
                T = Tn
            else:
                success = 0
                dT /= 2
                if debug: print 'Track step reduced to %.17f; condition = %s'%(
                        dT,
                        self.condition(prev_Z))
                if dT < 2.0**(-16):
                    raise ValueError('Track failed: step size limit reached.')
        return Zn
    
class Glunomial:
    """
    A product of powers of linear terms z_i or (1-z_i), as appears on
    the left side of a gluing equation.  These are Laurent monomials:
    powers may be negative.  Instantiate with one of the triples
    returned by Manifold.gluing_equations('rect').
    """
    def __init__(self, A, B, c):
        self.A, self.B, self.sign = array(A), array(B), float(c)

    def __repr__(self):
        apower = lambda n, p: 'z%d^%s'%(n,p) if p != 1 else 'z%s'%n
        bpower = lambda n, p: '(1-z%d)^%s'%(n,p) if p != 1 else '(1-z%s)'%n
        Apowers = [apower(n, a) for n, a in enumerate(self.A) if a != 0]
        Bpowers = [bpower(n, b) for n, b in enumerate(self.B) if b != 0]
        sign = '' if self.sign == 1.0 else '-'
        return sign + '*'.join(Apowers+Bpowers)

    def __call__(self, Z):
        W = 1 - Z
        try:
            return self.sign*prod(Z**self.A)*prod(W**self.B)
        except ValueError:
            print 'Glunomial evaluation crashed on %s'%self
            print 'A =', self.A
            print 'B =', self.B
            print 'c =', self.sign
            print 'Z =', Z
            raise ValueError
    def gradient(self, Z):
        W = 1 - Z
        return self.sign*prod(Z**self.A)*prod(W**self.B)*(self.A/Z - self.B/W)

class ShapeVector:
    """
    A vector of shape parameters, stored as a numpy.array, but
    with an approximate equality operator and distance operator ^.
    Instantiate with a sequence of complex numbers.  To access
    the underlying array, call as a function.
    """
    def __init__(self, values):
        self.array = array(values)

    def __str__(self):
        return self.array.__str__()

    def __getitem__(self, value):
        return self.array[value]

    def __len__(self):
        return len(self.array)
        
    def __eq__(self, other):
        return norm(self.array - other.array) < 1.0E-6

    def dist(self, other):
        return norm(self.array - other.array)

    def __call__(self):
        return self.array

    def __repr__(self):
        return repr(list(self))
    
    def is_degenerate(self):
        moduli = abs(self.array)
        return ( (moduli < 1.0E-6).any() or
                 (moduli < 1.0E-6).any() or
                 (moduli > 1.0E6).any()
                 )
    
class Fiber:
    """
    A fiber for the rational function [holonomy of the meridian]
    restricted to the curve defined by a PHCsystem.
    Manages a single PHCSytem with a complete solution set.  Can
    be initialized with a list of PHCsolutions.
    """
    def __init__(self, H_meridian, PHCsystem=None,
                 solutions=None, tolerance=1.0E-05):
        # Here the tolerance is used to determine which of the PHC solutions
        # are at infinity.
        self.H_meridian = H_meridian
        self.tolerance = tolerance
        if solutions:
            self.shapes = [ShapeVector(s) for s in solutions]
        self.system = PHCsystem
        if self.system:
            self.extract_info()
        
    def extract_info(self):
        N = self.system.num_variables()/2
        self.solutions = self.system.solution_list(tolerance=self.tolerance)
        # We only keep the "X" variables.
        self.shapes = [ShapeVector(S.point[:N]) for S in self.solutions]

    def __repr__(self):
        return 'Fiber(%s,\nsolutions=%s\n)'%(
            repr(self.H_meridian),
            repr([list(x) for x in self.shapes]).replace('],','],\n')
            )
    
    def __len__(self):
        return len(self.shapes)

    def __get_item(self, index):
        return self.shapes[index]
    
    def __eq__(self, other):
        """
        This ignores multiplicities.
        """
        for p in self.shapes:
            if p not in other.shapes:
                return False
        for p in other.shapes:
            if p not in self.shapes:
                return False
        return True
    
    def collision(self):
        for n, p in enumerate(self.shapes):
            for q in self.shapes[n+1:]:
                if p.dist(q) < 1.0E-10:
                    return True
        return False

    def is_finite(self):
        """
        Check if any cross-ratios are 0 or 1
        """
        for p in self.shapes:
            if p.is_degenerate():
                return False
        return True
            
    def details(self):
        # broken if not instantiated with a PHCsystem
        for n, s in enumerate(self.solutions):
            print 'solution #%s:'%n
            print s

    def residuals(self):
        # broken if not instantiated with a PHCsystem
        for n, s in enumerate(self.solutions):
            print n, s.res 

    def polish(self):
        # broken if not instantiated with a PHCsystem
        if self.system:
            self.system.polish()
            self.extract_info()

    def Tillmann_points(self):
        # broken if not instantiated with a PHCsystem
        if self.system is None:
            return []
        result = []
        for n, s in enumerate(self.solutions):
            if (s.t != 1.0 or self.shapes[n].is_degenerate()):
                result.append(n)
        return result

    def permutation(self, other):
        """
        return a list of pairs (m, n) where self.shapes[m] is
        closest to other.shapes[n].
        """
        result = []
        target = set(range(len(other.shapes)))
        for m, shape in enumerate(self.shapes):
            dist, n = min([(shape.dist(other.shapes[k]), k,)
                           for k in target])
            result.append( (m, n) )
            target.remove(n)
        return result

class PHCFibrator:
    """
    A factory for Fibers, computed by PHC or by a GluingSystem
    """
    def __init__(self, mfld, radius=1.02, saved_base_fiber=None, tolerance=1.0E-5):
        # Here the tolerance is used to determine when PHC solutions are
        # at infinity.
        self.manifold = mfld
        self.mfld_name = mfld.name()
        self.gluing_system = GluingSystem(mfld)
        self.radius = radius
        self.tolerance=tolerance
        self.num_tetrahedra = N = self.manifold.num_tetrahedra()
        variables = ( ['X%s'%n for n in range(N)] +
                      ['Y%s'%n for n in range(N)] )
        self.ring = PolyRing(variables + ['t'])
        self.equations = self.build_equations()
        self.equations += ['X%s + Y%s - 1'%(n,n) for n in range(N)] 
        self.psystem = ParametrizedSystem(
            self.ring,
            't',
            [PHCPoly(self.ring, e) for e in self.equations]
            )
        if saved_base_fiber and os.path.exists(saved_base_fiber):
            print 'Loading the starting fiber from %s'%saved_base_fiber
            datafile = open(saved_base_fiber)
            data = datafile.read()
            datafile.close()
            self.base_fiber = eval(data)
            self.basepoint = self.base_fiber.H_meridian
        else:
            self.basepoint = radius*exp(2*pi*1j*random())
            print 'Computing the starting fiber ... ',
            begin = time.time()
            self.base_system = self.psystem.start(self.basepoint, self.tolerance)
            print 'done. (%s seconds)'%(time.time() - begin)
            self.base_fiber = Fiber(self.basepoint,
                                    self.base_system)
            if saved_base_fiber:
                datafile = open(saved_base_fiber, 'w')
                datafile.write(repr(self.base_fiber))
                datafile.close()
                print 'Saved as %s'%saved_base_fiber

    def __len__(self):
        return len(self.base_fiber.solutions)
    
    def rect_to_PHC(self, eqn, rhs=None):
        A, B, c = eqn
        left = []
        if rhs is None:
            right = []
        elif isinstance(rhs, str):
            right = [rhs]
        else:
            right = [str(complex(rhs)).replace('j','*i')]
        for n, a in enumerate(A):
            if a > 0:
                left += ['X%s'%n]*a
            else:
                right += ['X%s'%n]*(-a)
        for n, b in enumerate(B):
            if b > 0:
                left += ['Y%s'%n]*b
            else:
                right += ['Y%s'%n]*(-b)
        if len(left) == 0:
            left = ['1']
        if len(right) == 0:
            right = ['1']
        op = ' - ' if c == 1 else ' + '
        return '*'.join(left) + op + '*'.join(right)

    def build_equations(self):
        if (self.manifold.num_cusps() != 1 or
            not self.manifold.is_orientable()):
            raise ValueError('Manifold must be orientable with one cusp.')
        eqns = self.manifold.gluing_equations('rect')
        meridian = eqns[-2]
        result = []
        for eqn in eqns[:-3]:
            result.append(self.rect_to_PHC(eqn))
        result.append(self.rect_to_PHC(meridian, rhs='t'))
        return result

    def PHC_to_complex(self, line):
       var, colon, real, imag = line.split()
       index=int(var[1:])
       if imag[0] == '-':
	   op = ''
       else:
	   op = '+'
       try:
	   return var[0], index, complex('%s%s%sj'%(real,op,imag))
       except ValueError:
	   print 'PHC parse error on %s+%sj'%(real,imag)

    def transport(self, start_fiber, target_holonomy, allow_collisions=False):
        """
        Use PHC to transport fibers.
        """
        target_system = self.psystem.transport(start_fiber.system,
                                               target_holonomy,
                                               allow_collisions)
        return Fiber(target_holonomy, target_system)

    def transport2(self, start_fiber, target_holonomy, allow_collisions=False,
                   debug=False):
        """
        Use a GluingSystem to transport fibers.
        """
        solutions = []
        dT = 1.0
        while True:
            if dT < 1.0/64:
                raise ValueError('Collision unavoidable. Try a different radius.')
            for shapevector in start_fiber.shapes:
                Zn = self.gluing_system.track(shapevector(),
                                              target_holonomy,
                                              dT=dT,
                                              debug=debug)
                solutions.append(Zn)
            result = Fiber(target_holonomy, solutions=solutions)
            if result.collision():
                dT *= 0.5
            else:
                break
        return result
    
class Holonomizer:
    """
    A family of fibers for the meridian holonomy map, lying
    above the points Rξ^m where ξ =e^(2πi/N). (N=128 by default.)
    The holonomizer can be tightened, to also contain fibers over
    points on the circle of radius T (= 1.0 by default).
    """
    def __init__(self, manifold, order=128, radius=1.02, saved_base_fiber=None):
        self.order = order
        self.radius = radius
        self.fibrator = PHCFibrator(manifold, radius, saved_base_fiber)
        self.manifold = manifold
        self.hp_manifold = self.manifold.high_precision()
        self.betti2 = [c % 2 for c in manifold.homology().coefficients].count(0)
        self.base_fiber = self.fibrator.base_fiber
        if not self.base_fiber.is_finite():
            raise RuntimeError, 'The starting fiber contains Tillmann points.'
        self.degree = len(self.base_fiber)
        print 'Degree is %s.'%self.degree
        # pre-initialize by just inserting an integer for each fiber
        # if the fiber construction fails, this can be detected by
        # isinstance(fiber, int)
        self.R_fibers = range(order)
        self.T_fibers = range(order)
        self.dim = manifold.num_tetrahedra()
        self.rhs = []
        eqns = manifold.gluing_equations('rect')
        self.glunomials = [Glunomial(A, B, c) for A, B, c in eqns[:-3]]
        self.rhs = [1.0]*(len(eqns) - 3)
        self.M_holo, self.L_holo = [Glunomial(A,B,c) for A,B,c in eqns[-2:]]
        self.glunomials.append(self.M_holo)
        start = time.time()
        self.track_satellite()
        print 'tracked in %s seconds.'%(time.time() - start)
        self.R_longitude_holos, self.R_longitude_evs = self.longidata(
            self.R_fibers)
            
    def __call__(self, Z):
        return array([F(Z) for F in self.glunomials])

    def track_satellite(self):
        """
        Construct the fibers over the circle of radius R.
        """
        print 'Tracking the satellite at radius %s ...'%self.radius
        arg = log(self.base_fiber.H_meridian).imag%(2*pi)
        R = self.radius
        Darg = 2*pi/self.order
        # The minus makes us consistent with the sign convention of numpy.fft
        self.R_circle = circle = [R*exp(-n*Darg*1j) for n in range(self.order)]
        self.base_index = base_index = (self.order - int((arg)/Darg))%self.order
        print ' %-5s\r'%base_index,
        # Move to the R-circle, if necessary.
        self.R_fibers[base_index] = self.fibrator.transport2(
                self.base_fiber, circle[base_index])
        for n in xrange(base_index+1, self.order):
            print ' %-5s\r'%n,
            sys.stdout.flush()
            self.R_fibers[n] = F = self.fibrator.transport2(
                self.R_fibers[n-1], circle[n])
            # self.R_fibers[n].polish()
            if not F.is_finite():
                print '**',
        for n in xrange(base_index-1, -1, -1):
            print ' %-5s\r'%n,
            sys.stdout.flush()
            self.R_fibers[n] = F = self.fibrator.transport2(
                self.R_fibers[n+1], circle[n])
            if not F.is_finite():
                print '**',
        print
        self.last_R_fiber = self.fibrator.transport2(self.R_fibers[-1],
                                                    self.R_fibers[0].H_meridian)
        print 'Polishing the end fibers ...'
        self.R_fibers[0].polish()
        self.last_R_fiber.polish()
        print 'Checking for completeness ... ',
        if not self.last_R_fiber == self.R_fibers[0]:
            print 'The end fibers did not agree!'
            print 'It might help to use a larger radius, or you might'
            print 'have been unlucky in your choice of base fiber.'
        else:
            print 'OK'

    def wiggle_track(self, R=1.02):
        """
        Construct the fibers over the unit circle, trying to avoid
        collisions by bumping out at singularities.
        """
        print 'Wiggle tracking'
        self.wiggle_fibers = [None for n in range(self.order)]
        bumps = []
        R = 1.0
        arg = log(self.base_fiber.H_meridian).imag%(2*pi)
        Darg = 2*pi/self.order
        # The minus makes us consistent with the sign convention of numpy.fft
        circle = [exp(-n*Darg*1j) for n in range(self.order)]
        base_index = (self.order - int(arg/Darg))%self.order
        print 'base index is %d'%base_index
        print ' %-5s\r'%base_index,
        self.wiggle_fibers[base_index] = self.fibrator.transport2(
                self.base_fiber, circle[base_index])
        for n in xrange(base_index+1, self.order):
            print ' %-5s\r'%n,
            sys.stdout.flush()
            F = self.wiggle_fibers[n-1]
            try:
                F = self.fibrator.transport2(F, circle[n])
            except ValueError:
                print '\nbumping at %d'%n
                bumps.append(n)
                F = self.fibrator.transport2(F, R*circle[n-1])
                F = self.fibrator.transport2(F, R*circle[n])
            # self.R_fibers[n].polish()
            if not F.is_finite():
                print '**',
            self.wiggle_fibers[n] = F
        print '\nreversing'
        F = self.wiggle_fibers[base_index]
        for n in xrange(base_index-1, -1, -1):
            print ' %-5s\r'%n,
            sys.stdout.flush()
            try:
                F = self.fibrator.transport2(F, circle[n])
            except ValueError:
                F = self.fibrator.transport2(F, 1.1*circle[n+1])
                F = self.fibrator.transport2(F, 1.1*circle[n])
            if not F.is_finite():
                print '**',
            self.wiggle_fibers[n] = F
        print
        self.last_wiggle_fiber = self.fibrator.transport2(
            self.wiggle_fibers[-1], self.wiggle_fibers[0].H_meridian)
        print 'Polishing the end fibers ...'
        self.wiggle_fibers[0].polish()
        self.last_wiggle_fiber.polish()
        print 'Checking for completeness ... ',
        if not self.last_wiggle_fiber == self.wiggle_fibers[0]:
            print 'The end fibers did not agree!'
            print 'It might help to use a larger radius, or you might'
            print 'have been unlucky in your choice of base fiber.'
        else:
            print 'OK'
        print 'tightening'
        for n in bumps:
            try:
                self.wiggle_fibers[n] = self.fibrator.transport2(
                    self.wiggle_fibers[n], circle[n])
            except ValueError:
                print 'failed to smooth bump at %s'%n

    def tighten(self, T=1.0):
        print 'Tightening the circle to radius %s ...'%T
        Darg = 2*pi/self.order
        self.T_circle = circle = [T*exp(-n*Darg*1j) for n in range(self.order)]
        for n in xrange(self.order):
            print ' %-5s\r'%n,
            sys.stdout.flush()
            try:
                self.T_fibers[n] = self.fibrator.transport2(self.R_fibers[n], circle[n])
            except ValueError:
                print 'Tighten failed at %s'%n
            #self.T_fibers[n] = self.fibrator.transport(self.R_fibers[n], circle[n])
            #self.T_fibers[n].system.polish()
        print '\nChecking for Tillmann points.'
        for n in xrange(self.order):
            try:
                t = self.T_fibers[n].Tillmann_points()
                if t:
                    print 'Tillmann points %s found in fiber %s.'%(t, n)
            except AttributeError: # If the fiber was not computed.
                print 'skipping %s'%n
        self.T_longitude_holos, self.T_longitude_evs = self.longidata(self.T_fibers)

    def longidata(self, fiber_list):
        """Compute the longitude holonomies and eigenvalues at each point in
        each fiber in the argument.  We allow the list to contain
        placeholders for failed computations.  A placeholder is any
        object which is not an instance of Fiber.  The fibers must all
        have degree == self.degree.

        The method produces two lists of lists.  Each outer list has
        length self.degree.  The inner lists contain pairs (n,z) where
        n is the index of the fiber in the input list and z is the
        associated holonomy or eigenvalue.  The integer can be used to
        reconstruct the corresponding meridian holonomy or eigenvalue
        in the case where the list of fibers is the elevation of a
        sampled circle.

        """
        print 'Computing longitude holonomies and eigenvalues.'
        # This crashes if there are bad fibers.
        longitude_holonomies = [
            [( n, self.L_holo(f.shapes[m]()) ) for n, f in enumerate(fiber_list)
             if isinstance(f, Fiber)]
            for m in xrange(self.degree)]
        # I tried choosing a random fiber here, and things broke badly.
        # find_shift would get the wrong shift.
#        if index is None:
#            index = randint(0,self.order - 1)
        index = 0 if isinstance(fiber_list[0], Fiber) else randint(0,self.order - 1)
        print 'starting index = %d'%index 
        longitude_traces = self.find_longitude_traces(fiber_list[index])
        longitude_eigenvalues = []
        for m, L in enumerate(longitude_holonomies):
            tr = longitude_traces[m]
            n, holo = L[index]
            e = sqrt(holo)
            # Choose the sign for the eigenvalue at the starting fiber
            E = [ (n,e) if abs(e + 1/e - tr) < abs(e + 1/e + tr) else (n,-e) ]
            # Avoid discontinuities caused by the branch cut used by sqrt
            for n, holo in L[index+1:]:
                e = sqrt(holo)
                E.append( (n,e) if abs(e - E[-1][1]) < abs(e + E[-1][1]) else (n,-e) )
            if index > 0:
                for n, holo in L[index-1::-1]:
                    e = sqrt(holo)
                    E.insert(0,(n,e) if abs(e - E[0][1]) < abs(e + E[0][1]) else (n,-e) )
            longitude_eigenvalues.append(E)
        return longitude_holonomies, longitude_eigenvalues

    def permutation(self, fiber_list):
        result = Permutation()
        start, end = fiber_list[0].shapes, fiber_list[-1].shapes
        for n, shape in enumerate(start):
            distances = array([shape.dist(end_shape) for end_shape in end])
            result[n] = argmin(distances)
        return result

    def compute_volumes(self, fiber_list):
        volumes = [ [] for n in range(self.degree) ]
        for fiber in fiber_list:
            for n, shape in enumerate(fiber.shapes):
                self.manifold.set_tetrahedra_shapes(shape(), fillings=[(0,0)])
                volumes[n].append(self.manifold.volume())
        return volumes
        
    def find_longitude_traces(self, fiber):
        # Sage complex numbers do not support attributes .real and .imag :^(((
        trace = lambda rep : complex(rep[0,0] + rep[1,1])
        traces = []
        for shape in fiber.shapes:
            sh = shape()
            self.hp_manifold.set_tetrahedra_shapes(sh, sh, [(0,0)])
            G = self.hp_manifold.fundamental_group()
            longitude = G.peripheral_curves()[0][1]
            relators = G.relators()
            generators = G.generators()
            M, N = len(relators), G.num_generators()
            A = matrix(zeros((M,N),'i'))
            L = zeros(N,'i')
            rhs = zeros(M,'i')
            for i in range(M):
                for j in range(N):
                    A[i,j] = (relators[i].count(generators[j]) +
                               relators[i].count(generators[j].upper()))%2
                    L[j] = (longitude.count(generators[j]) +
                               longitude.count(generators[j].upper()))%2
                rhs[i] = trace(G.SL2C(relators[i])).real < 0
            S = matrix(solve_mod2_system(A,rhs)).transpose()
            # Paranoia
            if max((A*S - matrix(rhs).transpose())%2) > 0:
                if self.betti2 == 1:
                    raise RuntimeError, "Mod 2 solver failed!"
            tr = trace(G.SL2C(longitude))
            if (L*S)%2 != 0:
                tr = -tr
            traces.append(tr)
        return traces

    def SL2C(self, word, shape):
        S = shape()
        self.hp_manifold.set_tetrahedra_shapes(S, None, [(0,0)])
        G = self.hp_manifold.fundamental_group()
        return G.SL2C(word)

    def O31(self, word, shape):
        S = shape()
        self.hp_manifold.set_tetrahedra_shapes(S, None, [(0,0)])
        G = self.hp_manifold.fundamental_group()
        return G.O31(word)

    def has_real_traces(self, shape):
        tolerance = 1.0E-10
        gens = self.hp_manifold.fundamental_group().generators()
        gen_mats = [self.SL2C(g, shape) for g in gens]
        for A in gen_mats:
            tr = complex(A[0,0] + A[1,1])
            if abs(tr.imag) > tolerance:
                return False
        mats = gen_mats[:]
        for i in range(1, len(gens) + 1):
            new_mats = []
            for A in gen_mats:
                for B in mats:
                    C = B*A
                    tr = complex(C[0,0] + C[1,1])
                    if abs(tr.imag) > tolerance:
                        return False
                    new_mats.append(C)
        return True
        
    def in_SU2(self, shape):
        tolerance = 1.0E-5
        gens = self.hp_manifold.fundamental_group().generators()
        # Check that all generators have real trace in [-2,2]
        for S in [self.SL2C(g, shape) for g in gens]:
            tr = complex(S[0,0] + S[1,1])
            if abs(tr.imag) > tolerance:
#                print 'trace not real'
                return False
            if abs(tr.real) > 2.0:
#                print 'trace not in [-2,2]'
                return False
        # Get O31 matrix generators
        o31matrices = [real_array(array(self.O31(g, shape))) for g in gens]
        # Take the first two
        A, B = o31matrices[:2]
        # find their axes
        M = matrix(zeros((4,4)))
        u, s, v = svd(A - eye(4))
        vt = transpose(v)
        M[:,[0,1]] = vt[:,[n for n in range(4) if abs(s[n]) < tolerance]]
        u, s, v = svd(B - eye(4))
        vt = transpose(v)
        M[:,[2,3]] = vt[:,[n for n in range(4) if abs(s[n]) < tolerance]]
        # Check if the axes cross, and find the fixed point (i.e. Minkwoski line)
        u, s, v = svd(M)
        vt = transpose(v)
        rel = vt[:,[n for n in range(4) if abs(s[n]) < tolerance]]
        if rel.shape != (4,1):
#            print 'linear algebra failure'
            return False
        # We have two descriptions -- average them
        rel[2] = -rel[2]
        rel[3] = -rel[3]
        fix = M*rel
        # check if the fixed line is in the light cone
        if abs(fix[0]) <= norm(fix[1:]):
#            print 'fixed line is not in the light cone'
            return False
        # check if all generators fix the same point
        for O in o31matrices:
            if norm(O*fix - fix) > tolerance:
#                print 'some generators do not share the fixed point.'
                return False
        return True
    
    def show_R_longitude_evs(self):
        R_plot = Plot([[x for n,x in track] for track in self.R_longitude_evs])

    def show_T_longitude_evs(self):
        T_plot = Plot([[x for n,x in track] for track in self.T_longitude_evs])

    def holo_permutation(self):
        return [self.R_fibers[0].shapes.index(p)
                for p in self.last_R_fiber.shapes]

    def holo_orbits(self):
        P = self.holo_permutation()
        Q = list(P)
        orbits = []
        while Q:
            first_one = this_one = Q.pop()
            orbit = []
            while True:
                orbit.append(this_one)
                this_one = P[this_one]
                if this_one == first_one:
                    break
                else:
                    Q.remove(this_one)
            orbits.append(orbit)
        return orbits

def solve_mod2_system(the_matrix,rhs):
    M,N = the_matrix.shape
    A = zeros((M,N+1),'i')
    A[:,:-1] = the_matrix
    A[:,-1] = rhs
    S = zeros(N,'i')
    P = []
    R = range(M)
    r = 0
    for j in range(N):
        i = r
        while i < M:
            if A[R[i]][j] != 0:
                break
            i += 1
        if i == M:
            continue
        if i > r:
            R[r], R[i] = R[i], R[r]
        P.insert(0,j)
        for i in range(r+1,M):
            if A[R[i]][j] == 1:
                A[R[i]] = A[R[i]]^A[R[r]]
        r += 1
    i = len(P)-1 
    for j in P:
        S[j] = (A[R[i]][N] - dot(A[R[i]][j+1:-1], S[j+1:]))%2
        i -= 1
    return S

class PEArc(list):
    """
    A list of pillowcase points lying on an arc of the PECharVariety.
    Subclassed here to allow additional attributes.
    """

class PECharVariety:
    def __init__(self, mfld, order=128, radius=1.02,
                 holonomizer=None, base_dir='base_fibers', hint_dir='hints'):
        if isinstance(mfld, Manifold):
            self.manifold = mfld
            self.manifold_name = mfld.name()
        else:
            self.manifold = Manifold(mfld)
            self.manifold_name = mfld
        saved_base_fiber = os.path.join(base_dir,
                                        self.manifold.name()+'.base')
        if holonomizer is None:
            self.holonomizer = Holonomizer(
                self.manifold,
                order=order,
                radius=radius,
                saved_base_fiber=saved_base_fiber)
            self.holonomizer.tighten()
        else:
            self.holonomizer = holonomizer
        self.order = self.holonomizer.order
        self.manifold = self.holonomizer.manifold

    # not used
    def build_components(self):
        H = self.holonomizer
        delta_M = -1.0/self.order
        M_args = 0.5 * ( arange(self.order, dtype=float64)*delta_M % 1.0 )
        # Construct the subarcs which lie on the unit torus
        arcs = []
        for track in self.holonomizer.T_longitude_evs:
            arc = []
            for n, ev in enumerate(track):
                if (0.99999 < abs(ev) < 1.00001):
                    arc.append( (log(ev).imag/(2*pi))%1.0 + M_args[m]*1j )
            arcs.append(arc)
        # Group the subarcs into components and add the cups and caps
        orbits = H.permutation(H.R_fibers).orbits()
        self.components = [ [arcs[n] for n in orbit] for orbit in orbits]
        for component in self.components:
            # add caps
            component.sort(key=lambda x : x[0].imag, reverse=True)
            n = len(component)
            while n > 0:
                if component[n][0].imag == component[n-1][0].imag:
                    component[n].insert(0, component[n-1][0])
                    n -= 1
                n -= 1
            # add cups
            component.sort(key=lambda x : x[-1].imag)
            n = len(component)
            while n > 0:
                if component[n][-1].imag == component[n-1].imag:
                    component[n].append(component[n-1][-1])
                    n -= 1
                n -= 1

    def build_arcs(self, su2_only=False):
        self.arcs = []
        self.arc_info = []
        H = self.holonomizer
        delta_M = -1.0/self.order
        M_args = 0.5 * ( arange(self.order, dtype=float64)*delta_M % 1.0 )
        for m, track in enumerate(self.holonomizer.T_longitude_evs):
            arc, info = PEArc(), []
            for n, ev in track:
                su2_ok = True
                if su2_only:
                    try:
                        su2_ok = H.in_SU2(H.T_fibers[n].shapes[m])
                    except:
                        #For now, just throw it in so we can look at it.
                        print 'Exception in SU_2 test.'
                        su2_ok = True
                if (0.99999 < abs(ev) < 1.00001) and su2_ok:
                    L = (log(ev).imag/(2*pi))%1.0
                    if len(arc)>2:  # don't do this near the corners.
                        last_L = arc[-1].real
                        if last_L > 0.8 and L < 0.2:   # L became > 1
                            length = 1.0 - last_L + L 
                            interp = ((1.0-last_L)*M_args[n] + L*M_args[n-1])/length
                            arc.append(PEPoint(1.0, interp, leave_gap=True)) 
                            arc.append(PEPoint(0.0, interp))
                        elif last_L < 0.2 and L > 0.8: # L became < 0
                            length = last_L + 1.0 - L 
                            interp = (last_L*M_args[n] + (1.0 - L)*M_args[n-1])/length
                            arc.append(PEPoint(0.0, interp, leave_gap=True))
                            arc.append(PEPoint(1.0, interp))
                    arc.append(PEPoint(L,M_args[n]))
                    info.append( (m,n) )
                else:
                    if len(arc) > 1:
                        m, n = arc.first_info = info[0]
                        arc.first_shape = self.holonomizer.T_fibers[n].shapes[m]
                        m, n = arc.last_info = info[-1]
                        arc.last_shape = self.holonomizer.T_fibers[n].shapes[m]
                        self.arcs.append(arc)
                        self.arc_info.append(info)
                    arc = PEArc()
                    info = []
            if arc:
                # Dumb repetition
                m, n = arc.first_info = info[0]
                arc.first_shape = self.holonomizer.T_fibers[n].shapes[m]
                m, n = arc.last_info = info[-1]
                arc.last_shape = self.holonomizer.T_fibers[n].shapes[m]
                self.arcs.append(arc)
                self.arc_info.append(info)
        self.curve_graph = curve_graph = Graph([], range(len(self.arcs)))
        self.add_extrema()
        # build the color dict
        self.colors = {}
        for n, component in enumerate(curve_graph.components()):
            for m in component:
                self.colors[m] = n
        # Clean up endpoints at the corners of the pillowcase.
        for arc in self.arcs:
            try:
                if abs(arc[1] - arc[0]) > 0.25:
                    if arc[0].imag > 0.45 and arc[1].imag  < 0.05 :
                        arc[0] = arc[0] - 0.5j
                    elif arc[0].imag < 0.05 and arc[1].imag > 0.45:
                        arc[0] = arc[0]+ 0.5j
                    if arc[0].real > 0.9 and arc[1].real < 0.1:
                        arc[0] = arc[0] - 1.0
                    elif arc[0].real < 0.1 and arc[1].real > 0.9:
                        arc[0] = arc[0] + 1.0
                if abs(arc[-1] - arc[-2]) > 0.25:
                    if abs(arc[-2].imag) < 0.05 and arc[-1].imag > 0.45:
                        arc[-1] = arc[-1] - 0.5j
                    elif arc[-2].imag > 0.45 and arc[-1].imag < 0.05:
                        arc[-1] = arc[-1] + 0.5j
                    if arc[-2].real < 0.1 and arc[-1].real > 0.9:
                        arc[-1] = arc[-1] - 1.0
                    elif arc[-2].real > 0.9 and arc[-1].real < 0.1:
                        arc[-1] = arc[-1] + 1.0
            except TypeError:
                pass

    def add_extrema(self):
        arcs = list(self.arcs)
        # caps
        arcs.sort(key=lambda x : x[0].imag, reverse=True)
        while arcs:
            arc = arcs.pop(0)
            level = [arc]
            while arcs and arc[0].imag == arcs[0][0].imag:
                    level.append(arcs.pop(0))
            while len(level) > 1:
                distances = array([level[0].first_shape.dist(a.first_shape)
                             for a in level[1:]])
                cap = [level.pop(0), level.pop(distances.argmin())]
                cap.sort(key=lambda a : a[0].real)
                left, right = cap
                if .05 < right[0].imag < .45:
                    join = True
                    if right[1].real > right[0].real and left[1].real < left[0].real:
                        right.insert(0, left[0])
                    elif right[1].real < right[0].real and left[1].real > left[0].real:
                        right.insert(0, PEPoint(1.0, left[0].imag))
                        left.insert(0, PEPoint(0.0, left[0].imag))
                    else:
                        join = False
                    if join:
                        self.curve_graph.add_edge(
                            self.arcs.index(left),
                            self.arcs.index(right))
        # cups
        arcs = list(self.arcs)
        arcs.sort(key=lambda x : x[-1].imag)
        while arcs:
            arc = arcs.pop(0)
            level = [arc]
            while arcs and arc[-1].imag == arcs[0][-1].imag:
                    level.append(arcs.pop(0))
            while len(level) > 1:
                distances = array([level[0].last_shape.dist(a.last_shape)
                             for a in level[1:]])
                cup = [level.pop(0), level.pop(distances.argmin())]
                cup.sort(key=lambda a : a[-1].real)
                left, right = cup
                if 0.05 < right[-1].imag < 0.45:
                    join = True
                    if right[-2].real > right[-1].real and left[-2].real < left[-1].real:
                        right.append(left[-1])
                    elif right[-2].real < right[-1].real and left[-2].real > left[-1].real:
                        right.append(PEPoint(1.0, left[-1].imag))
                        left.append(PEPoint(0.0, left[-1].imag))
                    else:
                        join = False
                    if join:
                        self.curve_graph.add_edge(
                            self.arcs.index(left),
                            self.arcs.index(right))                 
        return

    def XXXadd_extrema(self):
        arcs = list(self.arcs)
        # add caps
        arcs.sort(key=lambda x : x[0].imag, reverse=True)
        n = len(arcs) - 1
        while n > 0:
            if ( arcs[n][0].imag == arcs[n-1][0].imag
                 and abs(arcs[n][0].real - arcs[n-1][0].real) > 0.01 ):
                arcs[n].insert(0, arcs[n-1][0])
                self.curve_graph.add_edge(
                    self.arcs.index(arcs[n]),
                    self.arcs.index(arcs[n-1]))                 
                n -= 1
            n -= 1
        # add cups
        arcs.sort(key=lambda x : x[-1].imag)
        n = len(arcs) - 1
        while n > 0:
            if ( arcs[n][-1].imag == arcs[n-1][-1].imag > 0.01
                 and abs(arcs[n][-1].real - arcs[n-1][-1].real) > 0.01 ):
                arcs[n].append(arcs[n-1][-1])
                self.curve_graph.add_edge(
                    self.arcs.index(arcs[n]),
                    self.arcs.index(arcs[n-1]))                 
                n -= 1
            n -= 1
 
    def show(self, su2_only=False):
        self.build_arcs(su2_only)
        term = 'aqua' if sys.platform == 'darwin' else 'wxt'
        Plot(self.arcs,
             limits=((0.0, 1.0), (0.0, 0.5)),
             margins=(0,0),
             aspect='equal',
             title=self.manifold_name,
             linewidth=global_linewidth,
             colors = self.colors,
             extra_lines=[((0.5,0.5),(0.0,1.0))],
             extra_line_args={'color':'black', 'linewidth':0.75},
             commands="""
                    set terminal %s title "%s" size 1400,700
                    set xrange [0.0:1.0]
                    set yrange[0.0:0.5]
                    set xtics 0.5
                    set grid xtics nomxtics
                    set mxtics
                    """%(term, self.manifold_name),
             )
    
class Apoly:
    """
    The A-polynomial of a SnapPy manifold.  

    Constructor: Apoly(mfld, fft_size=128, gluing_form=False, denom=None, multi=False)
    <mfld>           is a manifold name recognized by SnapPy, or a Manifold instance.
    <gluing_form>    (True/False) indicates whether to find a "standard"
                     A-polynomial, or the gluing variety variant.
    <tight>           (True/False) if set, try to get the holonomizer to 
                     tighten to radius 1, and use those fibers.
    <fft_size>       must be at least twice the M-degree.  Try doubling this
                     if the coefficients seem to be wrapping.
    <denom>          Denominator for leading coefficient.  This should be
                     a string, representing a polynomial expression in H,
                     the meridian holonomy.  e.g. denom='(H-1)*(H-1)'
    <multi>          If true, multiple copies of lifts are not removed, so
                     multiplicities of factors of the polynomial are computed. 

  Methods:
    An Apoly object A is callable:  A(x,y) returns the value at (x,y).
    A.as_polynomial() returns a string suitable for input to a symbolic
                      algebra program.
    A.show_R_longitude_evs() uses gnuplot to graph the L-projections
                      of components of the inverse image of the satellite
                      circle of radius R in the M-plane.
    A.show_T_longitude_evs() uses gnuplot to graph the L-projections
                      of components of the inverse image of the tightened
                      circle of radius T in the M-plane.
    A.show_transforms() uses gnuplot to graph the coefficient polynomials
                      a_n(M), where A = a_0(M) + a_1(M)L + ...
    A.show_newton(text=False) shows the newton polygon with dots.  The text
                      flag shows the coefficients.
    A.boundary_slopes() prints the boundary slopes detected by the character
                      variety.
    A.save(basename=None, dir='polys', with_hint=True, twist=0)
                      Saves the polynomial in a .apoly or .gpoly text file for
                      input to a symbolic computation program.  The directory
                      can be overridden by specifying dir. Saves the parameters
                      in a .hint file unless with_hint==False.  Assumes that the
                      preferred longitude is LM^twist, where L,M are the SnapPea
                      meridian and longitued
    A.verify() runs various consistency checks on the polynomial.

    An Apoly object prints itself as a matrix of coefficients.
  """
    def __init__(self, mfld, fft_size=128, gluing_form=False, tight=False,
                 radius=1.02, denom=None, multi=False,
                 apoly_dir='apolys', gpoly_dir='gpolys',
                 base_dir='base_fibers', hint_dir='hints'):
        if isinstance(mfld, Manifold):
            self.manifold = mfld
            self.mfld_name = mfld.name()
        else:
            self.mfld_name = mfld
            self.manifold = Manifold(mfld)
        self.gluing_form = gluing_form
        self.tight = tight
        options = {'fft_size'    : fft_size,
                   'denom'       : denom,
                   'multi'       : multi,
                   'radius'      : radius,
                   'apoly_dir'   : apoly_dir,
                   'gpoly_dir'   : gpoly_dir,
                   'base_dir'    : base_dir,
                   'hint_dir'    : hint_dir}
#        if (fft_size, radius, denom, multi) == (128, None, None, False):
#            print "Checking for hints ...",
#            hintfile = os.path.join(self.hint_dir, mfld_name+'.hint')
#            if os.path.exists(hintfile):
#                print "yes!" 
#                exec(open(hintfile).read())
#                options.update(hint)
#            else:
#                print "nope."
        self.fft_size = N = options['fft_size']
        self.denom = options['denom']
        self.multi = options['multi']
        self.radius = options['radius']
        self.apoly_dir = options['apoly_dir']
        self.gpoly_dir = options['gpoly_dir']
        self.base_dir = options['base_dir']
        self.hint_dir = options['hint_dir']
        filename = self.manifold.name()+'.base'
        saved_base_fiber = os.path.join(self.base_dir, filename) 
        self.holonomizer = Holonomizer(
            self.manifold,
            order=self.fft_size,
            radius=self.radius,
            saved_base_fiber=saved_base_fiber)
        if self.tight:
            self.holonomizer.tighten(1.01)
            if self.gluing_form:
                vals = [array([x for n,x in track])
                        for track in self.holonomizer.R_longitude_holos]
            else:
                vals = [array([x for n,x in track])
                        for track in self.holonomizer.R_longitude_evs]
        else:
            if self.gluing_form:
                vals = [array([x for n,x in track])
                        for track in self.holonomizer.R_longitude_holos]
            else:
                vals = [array([x for n,x in track])
                        for track in self.holonomizer.R_longitude_evs]
        if multi == False:
            self.multiplicities, vals = self.demultiply(vals)
        self.sampled_coeffs = self.symmetric_funcs(vals)
        self.reduced_degree = len(self.sampled_coeffs)
        if self.denom:
            # denom is a polynomial in H = M^2 = holonomy of meridian.
            H = array(self.holonomizer.R_circle)
            exec('D = %s'%self.denom)
            self.raw_coeffs = array([ifft(x*D) for x in self.sampled_coeffs])
        else:
            self.raw_coeffs = array([ifft(x) for x in self.sampled_coeffs])
        # Renormalize, to account for R
        if N%2 == 0:
            renorm = self.radius**(-array(range(1+N/2)+range(1-N/2, 0)))
        else:
            renorm = self.radius**(-array(range(1+N/2)+range(-(N/2), 0)))
        self.normalized_coeffs = self.raw_coeffs*renorm
        self.int_coeffs = array([map(round, x.real)
                                 for x in self.normalized_coeffs])
        self.noise = self.normalized_coeffs.real - self.int_coeffs
        self.max_noise = [max(abs(x)) for x in self.noise]
        self.shift = self.old_find_shift()
        print 'Shift is %s'%self.shift
        if self.shift is None:
            print ('Coefficients may be wrapping.  '
                   'If so, a larger fft size would help.')
            return
        self.height = max([max(abs(x)) for x in self.int_coeffs])
        if self.height > float(2**52):
            print "Coefficients overflowed."
        C = self.int_coeffs.transpose()
        coefficient_array =  take(C, arange(len(C))-self.shift, axis=0)
        rows, cols = coefficient_array.shape
        while rows > 0:
            if max(abs(coefficient_array[rows-1])) > 0:
                break
            rows -= 1
        # This failed on 9_39 -- needs to remove empty rows
        self.coefficients = coefficient_array[:rows]
        print "Noise levels: "
        for level in self.max_noise:
            print level
        if max(self.max_noise) > 0.1:
            print 'Failed to find integer coefficients'
            return
        print 'Computing the Newton polygon.'
        power_scale = (1,1) if self.gluing_form else (1,2) 
        self.newton_polygon = NewtonPolygon(self.as_dict(), power_scale)
            
    def __call__(self, M, L):
        result = 0
        rows, cols = self.coefficients.shape
        for i in range(rows):
            Lresult = 0
            for j in range(cols):
                Lresult = Lresult*L + self.coefficients[-1-i][-1-j]
            result = result*M + Lresult
        return result
    
    def __repr__(self):
        digits = 2 + int(ceil(log(self.height)/log(10)))
        width = len(self.coefficients[0])
        format = '[' + ('%' + str(digits) + '.0f')*width + ']\n'
        result = ''
        for row in self.coefficients:
            result += format%tuple(row + 0.)
        return result

    def help(self):
        print self.__doc__

    def symmetric_funcs(self, roots):
        coeffs = [0, ones(roots[0].shape,'D')]
        for root in roots:
            for i in range(1, len(coeffs)):
                coeffs[-i] = -root*coeffs[-i] + coeffs[-1-i]
            coeffs.append(ones(roots[0].shape,'D'))
        return coeffs[1:]

    def demultiply(self, ev_list):
            multiplicities = []
            sdr = [] #system of distinct representatives
            multis = [1]*len(ev_list)
            for i in range(len(ev_list)):
                unique = True
                for j in range(i+1,len(ev_list)):
                    if max(abs(ev_list[i] - ev_list[j])) < 1.0E-6:
                        unique = False
                        multis[j] += multis[i]
                        break
                if unique:
                    sdr.append(i)
                    multiplicities.append((i, multis[i]))
            return multiplicities, [ev_list[i] for i in sdr]

    # this is the old version
    def XXfind_shift(self, raw_coeffs):
       rows, cols = raw_coeffs.shape
       N = self.fft_size
       shifts = [0]
       # Renormalize, since R != 1.
       if N%2 == 0:
           renorm = self.radius**(-array(range(1+N/2)+range(1-N/2, 0)))
       else:
           renorm = self.radius**(-array(range(1+N/2)+range(-(N/2), 0)))
       coeffs = raw_coeffs*renorm
       #start from the top and search for the last row above the middle
       #whose left-most non-zero entry is 1.
       for i in range(rows):
          for j in range(1, 1+ cols/2):
             if abs(abs(coeffs[i][-j]) - 1.) < .01:
                 shifts.append(j)
       print 'shifts: ', shifts
       return max(shifts), coeffs

    def old_find_shift(self):
       rows, cols = self.normalized_coeffs.shape
       shifts = [0]
       #start from the top and search for the last row above the middle
       #whose left-most non-zero entry is 1.
       for i in range(rows):
          for j in range(1, 1 + cols/2):
             if abs(abs(self.normalized_coeffs[i][-j]) - 1.) < .01:
                 shifts.append(j)
       print 'shifts (old method): ', shifts
       return max(shifts)

    # This is the new one, which fails on 9_3!!!
    def find_shift(self, cutoff=0.1):
        """
        Decide how many negative powers of M occur in the Laurent
        polynomial coefficients a_n(M).
        """
        N = self.fft_size
        coeffs = self.normalized_coeffs.transpose()
        # If the coefficients are large all the way to the middle, bail.
        if max(abs(coeffs[N/2])) > 0.5:
            return None
        # start from the bottom and go up until we hit a small row.
        for n in range(N-1, N/2, -1):
            if max(abs(coeffs[n])) < cutoff:
                return N-1-n 
        return None

# Should have a monomial class, and generate a list of monomials here not a string
    def monomials(self):
        rows, cols = self.coefficients.shape
        monomials = []
        for j in range(cols):
            for i in range(rows):
                if self.gluing_form:
                    m,n = 2*i, 2*j
                else:
                    m,n = 2*i, j
                a = int(self.coefficients[i][j])
                if a != 0:
                    if i > 0:
                        if j > 0:
                            monomial = '%d*(M^%d)*(L^%d)'%(a,m,n)
                        else:
                            monomial = '%d*(M^%d)'%(a,m)
                    else:
                        if j > 0:
                            monomial = '%d*(L^%d)'%(a,n)
                        else:
                            monomial = '%d'%a
                    monomials.append(monomial)
        return monomials

# Should use the list of monomials to generate the dict
    def as_dict(self):
        rows, cols = self.coefficients.shape
        result = {}
        for j in range(cols):
            for i in range(rows):
                if self.gluing_form:
                    m,n = 2*i, 2*j
                else:
                    m,n = 2*i, j
                coeff = int(self.coefficients[i][j])
                if coeff:
                    result[(m,n)] = coeff
        return result

    def break_line(self, line):
        marks = [0]
        start = 60
        while True:
            mark = line.find('+', start)
            if mark == -1:
                break
            marks.append(mark)
            start = mark+60
        lines = []
        for i in range(len(marks) - 1):
            lines.append(line[marks[i]:marks[i+1]])
        lines.append(line[marks[-1]:])
        return '\n    '.join(lines)
    
    def as_polynomial(self):
        polynomial = ('+'.join(self.monomials())).replace('+-','-')
        return polynomial

# could do this by sorting the monomials
    def as_Lpolynomial(self, name='A', twist=0):
        terms = []
        rows, cols = self.coefficients.shape
        #We are taking the true longitude to be L*M^twist.
        #So we change variables by L -> M^(-twist)*L.
        #Then renormalize so the minimal power of M is 0.
        minexp = 2*rows
        for j in range(cols):
            for i in range(rows):
                if self.coefficients[i][j]:
                    break
            minexp = min(2*i - j*twist, minexp)
        for j in range(cols):
            if self.gluing_form:
                n = 2*j
            else:
                n = j
            monomials = []
            for i in range(rows):
                m = 2*i
                a = int(self.coefficients[i][j])
                if a != 0:
                    if i > 0:
                        monomial = '%d*M^%d'%(a,m)
                    else:
                        monomial = '%d'%a
                    monomials.append(monomial.replace('^1 ',' '))
            if monomials:
                p = - n*twist - minexp
                if p:
                    P = '%d'%p
                    if p < 0:
                        P = '('+P+')'
                    if n > 0:
                        term = '+ (L^%d*M^%s)*('%(n,P) + ' + '.join(monomials) + ')'
                    else:
                        term = '(M^%s)*('%P + ' + '.join(monomials) + ')'
                else:
                    if n > 0:
                        term = '+ (L^%d)*('%n + ' + '.join(monomials) + ')'
                    else:
                        term = '(' + ' + '.join(monomials) + ')'
                term = self.break_line(term)
                terms.append(term.replace('+ -','- '))
        return name + ' :=\n' + '\n'.join(terms)

    def save(self, basename=None, dir=None, with_hint=True, twist=0):
        if dir == None:
            if self.gluing_form:
                poly_dir = self.gpoly_dir
                hint_dir = self.hint_dir
                ext = '.gpoly'
            else:
                poly_dir = self.apoly_dir
                hint_dir = self.hint_dir
                ext = '.apoly'
        for dir in (poly_dir, hint_dir):
            if not os.path.exists(dir):
                cwd = os.path.abspath(os.path.curdir)
                newdir = os.path.join(cwd,dir)
                response = raw_input("May I create a directory %s?(y/n)"%newdir)
                if response.lower()[0] != 'y':
                    sys.exit(0)
                os.mkdir(newdir)
        if basename == None:
            basename = self.mfld_name
        polyfile_name = os.path.join(poly_dir, basename + ext)
        hintfile_name = os.path.join(hint_dir, basename + '.hint')
        polyfile = open(polyfile_name,'w')
        if self.gluing_form:
            lhs = 'G_%s'%basename
        else:
            lhs = 'A_%s'%basename
        polyfile.write(self.as_Lpolynomial(name=lhs, twist=twist))
        polyfile.write(';\n')
        polyfile.close()
        if with_hint:
            hintfile = open(hintfile_name,'w')
            hintfile.write('hint={\n')
            hintfile.write('"radius" : %f,\n'%self.radius)
            hintfile.write('"fft_size" : %d,\n'%self.fft_size)
            if not self.multi:
                hintfile.write('"multi" : %s,\n'%self.multi)
            if self.denom:
                hintfile.write('"denom" : "%s",\n'%self.denom)
            hintfile.write('}\n')
            hintfile.close()
            
    def boundary_slopes(self):
        print self.newton_polygon.lower_slopes
        
    def show_R_longitude_evs(self):
        self.holonomizer.show_R_longitude_evs()

    def show_T_longitude_evs(self):
        self.holonomizer.show_T_longitude_evs()

    def show_coefficients(self):
        plot = Plot(self.normalized_coeffs.real)

    def show_noise(self):
        plot = Plot(self.noise)

    def show_imag_noise(self):
        plot = Plot(self.normalized_coeffs.imag)

    def show_newton(self, text=False):
        V = PolyViewer(self.newton_polygon, title=self.mfld_name)
        V.show_sides()
        if text:
            V.show_text()
        else:
            V.show_dots()

    def show_R_volumes(self):
        H = self.holonomizer
        Plot(H.compute_volumes(H.R_fibers))

    def show_T_volumes(self):
        H = self.holonomizer
        Plot(H.compute_volumes(H.T_fibers))

    def tighten(self):
        self.holonomizer.tighten()

    def verify(self):
        result = True
        sign = None
        print 'Checking noise level ...',
        print max(self.max_noise)
        if max(self.max_noise) > .3:
            result = False
            print 'Failed'
        print 'Checking for reciprocal symmetry ... ',
        maxgap = 0
        if max(abs(self.coefficients[0] - self.coefficients[-1][-1::-1]))==0:
            sign = -1.0
        elif max(abs(self.coefficients[0] + self.coefficients[-1][-1::-1]))==0:
            sign = 1.0
        else:
            print 'Failed!'
            result = False
        if sign:
            for i in range(len(self.coefficients)):
                maxgap = max(abs(self.coefficients[i] +
                              sign*self.coefficients[-i-1][-1::-1]))
                if maxgap > 0:
                    print 'Failed! gap = %d'%maxgap
                    result = False
        if result:
            print 'Passed!'
        return result

# This won't work with bad fibers.    
#    def tighten(self, T=1.0):
#      self.holonomizer.tighten(T)
#      roots = [array(x) for x in self.holonomizer.T_longitude_evs]
#      self.T_sampled_coeffs = self.symmetric_funcs(roots)
#      self.T_raw_coeffs = array([ifft(x) for x in self.T_sampled_coeffs])

class Slope:
      def __init__(self, xy, power_scale=(1,1)):
            x, y = xy
            x, y = x*power_scale[0], y*power_scale[1]
            if x == 0:
                  if y == 0:
                        raise ValueError('gcd(0,0) is undefined.')
                  else:
                        gcd = abs(y)
            else:
                  x0 = abs(x)
                  y0 = abs(y)
                  while y0 != 0:
                        r = x0%y0
                        x0 = y0
                        y0 = r
                  gcd = x0
            if x < 0:
                  x, y = -x, -y
            self.x = x/gcd
            self.y = y/gcd
      def __cmp__(self, other):
            return int.__cmp__(self.y*other.x - other.y*self.x, 0)
      def __repr__(self):
            return '%d/%d'%(self.y, self.x)

class NewtonPolygon:
      def __init__(self, coeff_dict, power_scale=(1,1)):
          # Clean up weird sage tuples
          self.coeff_dict = {}
          for (degree, coefficient) in coeff_dict.items():
              self.coeff_dict[tuple(degree)] = coefficient
          # The X-power is the y-coordinate!
          self.support = [(x[1], x[0]) for x in coeff_dict.keys()]
          self.support.sort()
          self.lower_slopes = []
          self.upper_slopes = []
          self.lower_vertices = []
          self.upper_vertices = []
          self.newton_sides = {}
          self.find_vertices()

      def slope(self, v, w):
          return Slope((w[0]-v[0], w[1]-v[1]))

      def find_vertices(self):
          last = self.support[0]
          T = []
          B = [last]
          for e in self.support[1:]:
              if e[0] != last[0]:
                  T.append(last)
                  B.append(e)
              last = e
          T.append(last)
          if T[0] != B[0]:
              self.lower_slopes.append(Slope((0,1)))
              self.upper_vertices.append(T[0])
          n = 0
          while n < len(B) - 1:
              self.lower_vertices.append(B[n])
              slopes = [(self.slope(B[n], B[k]), -k) for k in range(n+1,len(B))]
              slope, m = min(slopes)
              self.lower_slopes.append(slope)
              if slope.y < 0:
                  newton_side = [B[n]]
                  for s, j in slopes:
                      if s == slope and -j <= -m:
                          newton_side.append(B[-j])
                  self.newton_sides[(slope.x,slope.y)] = newton_side
              n = -m
          n = 0
          while n < len(T) - 1:
              slope, m = max([(self.slope(T[n], T[k]), k) for k in range(n+1,len(T))])
              self.upper_slopes.append(slope)
              self.upper_vertices.append(T[m])
              n = m
          if T[-1] != B[-1]:
              self.upper_slopes.append(Slope((0,1)))
              self.lower_vertices.append(B[-1])

      def side_dicts(self):
          result = {}
          for slope, side in self.newton_sides.items():
              side_dict = {}
              for i, j in side:
                  side_dict[(j,i)] = self.coeff_dict[(j,i)]
              result[slope] = side_dict
          return result
    
      def puiseux_expansion(self):
          result = []
          for slope, side_dict in self.side_dicts().items():
              P = sage_2poly(dict)
              m, n = slope
              result.append(P(t**n,t**m))



class PolyViewer:
      def __init__(self, newton_poly, title=None, scale=None, margin=50):
          self.NP = newton_poly
          self.columns = 1 + self.NP.support[-1][0]
          self.rows = 1 + max([d[1] for d in self.NP.support])
          if scale == None:
                scale = 600/max(self.rows, self.columns)
          self.scale = scale
          self.margin = margin
          self.width = (self.columns - 1)*self.scale + 2*self.margin
          self.height = (self.rows - 1)*self.scale + 2*self.margin
          self.window = Tkinter.Tk()
          if title:
              self.window.title(title)
          self.window.wm_geometry('+400+20')
          self.canvas = Tkinter.Canvas(self.window,
                                  bg='white',
                                  height=self.height,
                                  width=self.width)
          self.canvas.pack()
          self.font = ('Helvetica','18','bold')
          self.dots=[]
          self.text=[]
          self.sides=[]

          self.grid = (
              [ self.canvas.create_line(
              0, self.height - self.margin - i*scale,
              self.width, self.height - self.margin - i*scale,
              fill=self.gridfill(i))
                        for i in range(self.rows)] +
              [ self.canvas.create_line(
              self.margin + i*scale, 0,
              self.margin + i*scale, self.height,
              fill=self.gridfill(i))
                        for i in range(self.columns)])
#          self.window.mainloop()

      def write_psfile(self, filename):
            self.canvas.postscript(file=filename)
            
      def gridfill(self, i):
          if i:
              return '#f0f0f0'
          else:
              return '#d0d0d0'
          
      def point(self, pair):
          i,j = pair
          return (self.margin+i*self.scale,
                  self.height - self.margin - j*self.scale)
      
      def show_dots(self):
          r = 1 + self.scale/20
          for i, j in self.NP.support:
              x,y = self.point((i,j))
              self.dots.append(self.canvas.create_oval(
                  x-r, y-r, x+r, y+r, fill='black'))

      def erase_dots(self):
          for dot in self.dots:
              self.canvas.delete(dot)
          self.dots = []

      def show_text(self):
          r = 2 + self.scale/20
          for i, j in self.NP.support:
              x,y = self.point((i,j))
              self.sides.append(self.canvas.create_oval(
                  x-r, y-r, x+r, y+r, fill='black'))
              self.text.append(self.canvas.create_text(
                  2*r+x,-2*r+y,
                  text=str(self.NP.coeff_dict[(j,i)]),
                  font=self.font,
                  anchor='c'))
              
      def erase_text(self):
            for coeff in self.text:
                  self.canvas.delete(coeff)
            self.text=[]

      def show_sides(self):
            r = 2 + self.scale/20
            first = self.NP.lower_vertices[0]
            x1, y1 = self.point(first)
            upper = list(self.NP.upper_vertices)
            upper.reverse()
            vertices = self.NP.lower_vertices + upper + [first]
            for vertex in vertices:
                  self.sides.append(self.canvas.create_oval(
                        x1-r, y1-r, x1+r, y1+r, fill='red'))
                  x2, y2 = self.point(vertex)
                  self.sides.append(self.canvas.create_line(
                        x1, y1, x2, y2,
                        fill='red'))
                  x1, y1 = x2, y2

      def erase_sides(self):
            for object in self.sides:
                  self.canvas.delete(object)
            self.sides=[]

class Permutation(dict):
    def orbits(self):
        points = set(self.keys())
        orbits = []
        while points:
            first = n = points.pop()
            orbit = [first]
            while True:
                n = self[n]
                if n == first:
                    orbits.append(orbit)
                    break
                else:
                    points.remove(n)
                    orbit.append(n)
        return orbits

winding = lambda x : (sum(log(x[1:]/x[:-1]).imag) + log(x[0]/x[-1]).imag)/(-2*pi)
if got_sage:
    Apoly.sage_poly = lambda self : sage_2poly(self.as_dict(), ring=QQ['M','L'])
#    PolyRelation.sage_poly = lambda self : sage_poly(self.as_dict())

#M = Manifold('4_1')
#F = Fiber((-0.991020658402+0.133708842719j),
#          [
#           ShapeVector(array([
#            6.18394729421744E-01+5.14863122901458E-02j,
#            6.18394729421744E-01-5.14863122901458E-02j], dtype=DTYPE)),
#           ShapeVector(array([
#            -1.57365927858202E+00+3.47238981119960E-01j,
#            -1.57365927858202E+00-3.47238981119960E-01j], dtype=DTYPE))
#          ])
#begin = time.time()
#B = Holonomizer(M, F)
#print time.time()-begin
#Z = array(M.tetrahedra_shapes('rect'))
#print B.run_newton(Z, 1j)
