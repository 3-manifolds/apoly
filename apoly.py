from numpy import array, matrix, ndarray, dot, prod, diag, transpose, zeros, ones
from numpy import log, exp, pi, sqrt, ceil
from numpy import dtype, complex128, float64, take, arange, where
from numpy.linalg import svd, norm
from numpy.fft import ifft
try:
    from numpy import complex256 as big_complex
    from numpy import float128 as big_float
except ImportError:
    from numpy import complex192 as big_complex
    from numpy import float96 as big_float 
from phc import *
from snappy import *
from random import random
from subprocess import Popen, PIPE
import time, sys, os, Tkinter


DTYPE = dtype('c16')

class Glunomial:
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
        return self.sign*prod(Z**self.A)*prod((1 - Z)**self.B)

    def gradient(self, Z):
        return self(Z)*(self.A/Z - self.B/(1 - Z))

class Point:
    """
    A numpy array of complex numbers, with its own approximate
    equality operator.  Instantiate with a sequence of complex
    numbers.
    """
    def __init__(self, Z):
        self.Z = array(Z)
        
    def __repr__(self):
        return str(self.Z)

    def __eq__(self, other):
        return norm(self.Z - other.Z) < 1.0E-10

    def __xor__(self, other):
        return norm(self.Z - other.Z)

    def is_degenerate(self):
        return (min(abs(self.Z)) < 1.0E-6
                or min(abs(self.Z - 1.0)) < 1.0E-6
                or max(abs(self.Z)) > 1.0E6)
    
class Fiber:
    """
    A fiber for the rational function [holonomy of the meridian].
    Manages a single PHCSytem with a complete solution set.
    """
    def __init__(self, H_meridian, system, tolerance=1.0E-06):
        self.H_meridian = H_meridian
        self.system = system
        self.tolerance = tolerance
        self.extract_info()
        
    def extract_info(self):
        N = self.system.num_variables()/2
        self.solutions = self.system.solution_list()
        # We only keep the "X" variables.
        self.points = [Point(S.point[:N]) for S in self.solutions]

    def __len__(self):
        return len(self.solutions)

    def __get_item(self, index):
        return self.points[index]
    
    def __eq__(self, other):
        """
        This ignores multiplicities.
        """
        for p in self.points:
            if p not in other.points:
                return False
        for p in other.points:
            if p not in self.points:
                return False
        return True
    
    def collision(self):
        for n, p in enumerate(self.points):
            for q in self.points[n+1:]:
                if norm(p.Z - q.Z) < 1.0E-10:
                    return True
        return False

    def is_finite(self):
        """
        Check if any cross-ratios are 0 or 1
        """
        for p in self.points:
            if not p.is_degenerate:
                return False
        return True
            
    def details(self):
        for n, s in enumerate(self.solutions):
            print 'solution #%s:'%n
            print s

    def residuals(self):
        for n, s in enumerate(self.solutions):
            print n, s.res 

    def polish(self):
        self.system.polish()
        self.extract_info()

    def Tillmann_points(self):
        result = []
        for n, s in enumerate(self.solutions):
            if (s.t != 1.0 or self.points[n].is_degenerate()):
                result.append(n)
        return result
    
class PHCFibrator:
    """
    A factory for Fibers, computed by PHC.
    """
    def __init__(self, mfld_name, radius=1.02):
        self.mfld_name = mfld_name
        self.radius = radius
        self.manifold = Manifold(self.mfld_name)
        self.num_tetrahedra = N = self.manifold.num_tetrahedra()
        variables = ( ['X%s'%n for n in range(N)] +
                      ['Y%s'%n for n in range(N)] )
        self.ring = PolyRing(variables + ['t'])
        self.basepoint = radius*exp(2*pi*1j*random())
        self.equations = self.build_equations()
        self.equations += ['X%s + Y%s - 1'%(n,n) for n in range(N)] 
        self.psystem = ParametrizedSystem(
            self.ring,
            't',
            [PHCPoly(self.ring, e) for e in self.equations]
            )
        print 'Computing the starting fiber ... ',
        begin = time.time()
        self.base_system = self.psystem.start(self.basepoint, tolerance=1.0E-06)
        print 'done. (%s seconds)'%(time.time() - begin)
        self.base_fiber = Fiber(self.basepoint,
                                 self.base_system)

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
            raise ValueError, 'Manifold must be orientable with one cusp.'
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
        target_system = self.psystem.transport(start_fiber.system,
                                               target_holonomy,
                                               allow_collisions)
        return Fiber(target_holonomy, target_system)
    
class Holonomizer:
    """
    A family of fibers for the meridian holonomy map, lying
    above the Nth roots of unity on the unit circle.  (N=128 by default.)
    """
    def __init__(self, manifold_name, order=128, radius=1.02):
        self.order = order
        self.radius = radius
        self.fibrator = PHCFibrator(manifold_name, radius=radius)
        self.manifold = manifold = self.fibrator.manifold
        self.base_fiber = self.fibrator.base_fiber
        if not self.base_fiber.is_finite():
            raise RuntimeError, 'The starting fiber contains ideal points.'
        self.degree = len(self.base_fiber)
        self.dimension = manifold.num_tetrahedra()
        # pre-initialize
        self.R_fibers = range(order)
        self.T_fibers = range(order)
        self.dim = manifold.num_tetrahedra()
        self.rhs = []
        eqns = manifold.gluing_equations('rect')
        self.glunomials = [Glunomial(A, B, c) for A, B, c in eqns[:-3]]
        self.rhs = [1.0]*(len(eqns) - 3)
        self.M_holo, self.L_holo = [Glunomial(A,B,c) for A,B,c in eqns[-2:]]
        self.glunomials.append(self.M_holo)
        self.track_satellite()
        try:
            self.R_longitude_holos, self.R_longitude_evs = self.longidata(self.R_fibers)
        except:
            print 'Failed'
            
    def __call__(self, Z):
        return array([F(Z) for F in self.glunomials])

    def track_satellite(self): 
        print 'Tracking the satellite at radius %s ...'%self.radius
        arg = log(self.base_fiber.H_meridian).imag%(2*pi)
        R = self.radius
        Darg = 2*pi/self.order
        # The minus makes us consistent with the sign convention of numpy.fft
        self.R_circle = circle = [R*exp(-n*Darg*1j) for n in range(self.order)]
        self.base_index = base_index = (self.order - int((arg)/Darg))%self.order
        print base_index,
        self.R_fibers[base_index] = self.fibrator.transport(
                self.base_fiber, circle[base_index])
        for n in xrange(base_index+1, self.order):
            print n,
            sys.stdout.flush()
            self.R_fibers[n] = F = self.fibrator.transport(
                self.R_fibers[n-1], circle[n])
            self.R_fibers[n].system.polish()
            if not F.is_finite():
                print '**',
        for n in xrange(base_index-1, -1, -1):
            print n,
            sys.stdout.flush()
            self.R_fibers[n] = F = self.fibrator.transport(
                self.R_fibers[n+1], circle[n])
            if not F.is_finite():
                print '**',
        print
        self.last_R_fiber = self.fibrator.transport(self.R_fibers[-1],
                                                    self.R_fibers[0].H_meridian)
        print 'Polishing the end fibers ...'
        self.R_fibers[0].polish()
        self.last_R_fiber.polish()
        print 'Checking for completeness ... ',
        if not self.last_R_fiber == self.R_fibers[0]:
            print 'The end fibers did not agree!'
            print 'It might help to use a larger radius, or you might'
            print 'have been unlucky in your choise of base fiber.'
        else:
            print 'OK'

    def tighten(self, T=1.0):
        print 'Tightening the circle to radius %s ...'%T
        Darg = 2*pi/self.order
        self.T_circle = circle = [T*exp(-n*Darg*1j) for n in range(self.order)]
        for n in xrange(self.order):
            print n,
            sys.stdout.flush()
            self.T_fibers[n] = self.fibrator.transport(self.R_fibers[n], circle[n])
            self.T_fibers[n].system.polish()
        print '\nChecking for Tillmann points.'
        for n in xrange(self.order):
            t = self.T_fibers[n].Tillmann_points()
            if t:
                print 'Tillmann points %s found in fiber %s.'%(t, n)
        self.T_longitude_holos, self.T_longitude_evs = self.longidata(self.T_fibers)

    def longidata(self, fiber_list):
        print 'Computing longitude holonomies and eigenvalues.'
        longitude_holonomies = [
            [self.L_holo(f.points[n].Z) for f in fiber_list]
            for n in xrange(self.degree)]
        longitude_traces = self.find_longitude_traces(fiber_list[0])
        longitude_eigenvalues = []
        for n, L in enumerate(longitude_holonomies):
            tr = longitude_traces[n]
            E = []
            e = sqrt(L[0])
            E.append(e if abs(e + 1/e - tr) < abs(e + 1/e + tr) else -e ) 
            for m, H in enumerate(L[1:],1):
                e = sqrt(H)
                E.append(e if abs(e - E[m-1]) < abs(e + E[m-1]) else -e )
            longitude_eigenvalues.append(E)
        return longitude_holonomies, longitude_eigenvalues
    
    def find_longitude_traces(self, fiber):
        trace = lambda rep : rep[0,0] + rep[1,1]
        traces = []
        for point in fiber.points:
            #  I had to move the dehn_fill((0,0)) inside the loop to get this to
            #  work correctly when there is a denominator, e.g. for 8_17.
            #  The values of the longitude traces were coming out wrong.
            #  This was not needed with SnapPeaPython --- SnapPy bug???
            #  Why is it needed at all?
            self.manifold.dehn_fill((0,0))
            self.manifold.set_tetrahedra_shapes(point.Z, fillings=[(0,0)])
            G = self.manifold.fundamental_group()
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
                raise RuntimeError, "Mod 2 solver failed!"
            tr = trace(G.SL2C(longitude))
            if int((L*S)%2):
                tr = -tr
            traces.append(tr)
        return traces

    def jacobian(self, Z):
        return matrix([E.gradient(Z) for E in self.glunomials])

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

class SU2CharVariety:
    def __init__(self, manifold_name, order=128, radius=1.02, holonomizer=None):
        self.manifold_name = manifold_name
        if holonomizer is None:
            self.holonomizer = Holonomizer(manifold_name, order=order, radius=radius)
            self.holonomizer.tighten()
        else:
            self.holonomizer = holonomizer
        self.order = self.holonomizer.order
        self.manifold = self.holonomizer.manifold
        self.build_arcs()

    def build_arcs(self):
        self.arcs = []
        M_args = -arange(self.order, dtype=float64)/self.order
        for track in self.holonomizer.T_longitude_evs:
            arc = []
            for n, ev in enumerate(track):
                if 0.9999 < abs(ev) < 1.0001:
                    L = (1 + (log(ev).imag)/pi)%2
                    if n > 0 and abs(L - lastL) > 1 and len(arc) > 0:
                        arc.append(None)
                    arc.append( L + 1j*(M_args[n]%1) )
                    lastL = L
                else:
                    if len(arc) > 1:
                        self.arcs.append(arc)
                    arc = []
            if arc:
                self.arcs.append(arc)
        # Clean up endpoints at the corners of the pillowcase.
        for arc in self.arcs:
            try:
                if abs(arc[1] - arc[0]) > 0.8:
                    if abs(arc[0].imag) < .001:
                        arc[0] = arc[0] + 1j
                    elif abs(arc[0].imag - 1) < .001:
                        arc[0] = arc[0] - 1j
                if abs(arc[-1] - arc[-2]) > 0.8:
                    if abs(arc[-1].imag) < .001:
                        arc[-1] = arc[-1] + 1j
                    elif abs(arc[-1].imag - 1) < .001:
                        arc[-1] = arc[1] - 1j
            except TypeError:
                pass
                        
    def show(self):
        Plot(self.arcs, commands="""
                    set terminal aqua title "%s" size 1000 500
                    set for [i = 1:10] style line i lw 2
                    set xrange [0:2]
                    set yrange[0:1]
                    """%self.manifold_name)
    
class PolyRelation:
    """
    An integral polynomial relation satisfied by two coordinate
    functions M and L in the function field of a curve defined over Q.
    We view M as the "base element".  Then the polynomial relation is
    just a minimal polynomial for L over Q(M).  (In practice, for
    function fields of curves of characters, M is the holonomy of the
    meridian.)
    """
    def __init__(self, Lvalues, Mvalues, radius=1.02, fft_size=128,
                 denom=None, multi=False, gluing_form=False):
        self.radius = radius
        self.fft_size = fft_size
        self.denom = denom
        self.multi = multi
        self.gluing_form=gluing_form
        Larrays = [array(x) for x in Lvalues]
        if multi == False:
            self.multiplicities, Larrays = self.demultiply(Larrays)
        self.sampled_coeffs = self.symmetric_funcs(Larrays)
        if denom:
            M = array(Mvalues)
            exec('D = %s'%self.denom) #denom is an expression for a poly in M
            self.raw_coeffs = array([ifft(x*D) for x in self.sampled_coeffs])
        else:
            self.raw_coeffs = array([ifft(x) for x in self.sampled_coeffs])
        self.shift = self.find_shift(self.raw_coeffs)
        if self.shift is None:
            print 'Coefficients seem to be wrapping.  A larger fft size might help.'
            return
        renorm = self.radius**(-array(range(self.fft_size - self.shift)
                                      + range(-self.shift, 0)))
        self.float_coeffs = renorm*self.raw_coeffs
        self.height = max([max(abs(x.real)) for x in self.float_coeffs])
        # This has problems.
        if self.height > float(2**52):
            print "Coefficients overflowed."
        else:
            self.height = round(self.height)
        self.int_coeffs = array([map(round, x.real) for x in self.float_coeffs])
        C = self.int_coeffs.transpose()
        coefficient_array =  take(C, arange(len(C))-self.shift, axis=0)
        rows, cols = coefficient_array.shape
        while rows:
            if max(abs(coefficient_array[rows-1])) > 0:
                break
            rows -= 1
        self.coefficients = coefficient_array[:rows]
        self.newton_polygon = NewtonPolygon(self.coefficients, self.gluing_form)
        self.noise = [max(abs(self.float_coeffs[i] - self.int_coeffs[i])) for
                      i in range(len(self.float_coeffs))]
        print "Noise levels: "
        for level in self.noise:
            print level
            
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
            sdr = []
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

    def find_shift(self, raw_coeffs):
       rows, cols = raw_coeffs.shape
       N = self.fft_size
       shifts = [0]
       renorm = self.radius**(-array(range(1+N/2)+range(1-N/2, 0)))
       coeffs = raw_coeffs*renorm
       for i in range(rows):
          for j in range(1, 1+ cols/2):
             if abs(abs(coeffs[i][-j]) - 1.) < .001:
                 shifts.append(j)
       print 'shifts: ', shifts
       return max(shifts)


    def Xfind_shift(self, raw_coeffs, tolerance= 0.001):
       N = self.fft_size
       if N%2 == 0:
           renorm = self.radius**(-array(range(1+N/2)+range(1-N/2, 0)))
       else:
           renorm = self.radius**(-array(range(1+N/2)+range(-(N/2), 0)))
       coeffs = raw_coeffs*renorm
       maxes = [max(abs(row.real)) for row in coeffs.transpose()]
       for n in range(N):
           if maxes[-n] < tolerance:
               return n-1
    
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

    def as_polynomial(self):
        polynomial = ('+'.join(self.monomials())).replace('+-','-')
        return polynomial

    def show_newton(self, text=False):
        # The gluing_form should be renamed.  It indicates that the
        # polynomial is using the holonomies of the longitude, rather
        # than the eigenvalues.  In general, it would not make sense
        # to take square roots.
        V = Polyview(self.coefficients, self.gluing_form)
        V.show_sides()
        if text:
            V.show_text()
        else:
            V.show_dots()

#This should be derived from PolyRelation
class Apoly:
    """
    The A-polynomial of a SnapPy manifold.  

    Constructor: Apoly(mfld_name, fft_size=128, gluing_form=False, denom=None, multi=False)
    <mfld_name>      is a manifold name recognized by SnapPy.
    <gluing_form>    (True/False) indicates whether to find a "standard"
                     A-polynomial, or the gluing variety variant.
    <fft_size>       must be at least twice the M-degree.  Try doubling this
                     if the coefficients seem to be wrapping.
    <denom>          Denominator for leading coefficient.  This should be
                     a string, representing a polynomial expression in M.
    <multi>          If true, multiple copies of lifts are not removed, so
                     multiplicities of factors of the polynomial are computed. 

  Methods:
    An Apoly object A is callable:  A(x,y) returns the value at (x,y).
    A.as_polynomial() returns a string suitable for input to a symbolic
                      algebra program.
    A.show_lifts() uses gnuplot to graph the L-projections of components of
                   the inverse image of a circle in the M-plane.  The circle
                   is centered at the origin and the radius is as specified.
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
    def __init__(self, mfld_name, fft_size=128, gluing_form=False,
                 radius=1.02, denom=None, multi=False):
        self.mfld_name = mfld_name
        self.gluing_form = gluing_form
        options = {'fft_size'    : fft_size,
                   'denom'       : denom,
                   'multi'       : multi,
                   'radius'      : radius}
#        if (fft_size, radius, denom, multi) == (128, None, None, False):
#            print "Checking for hints ...",
#            hintfile = os.path.join(self.hint_dir, mfld_name+'.hint')
#            if os.path.exists(hintfile):
#                print "yes!" 
#                exec(open(hintfile).read())
#                options.update(hint)
#            else:
#                print "nope."
        self.fft_size = options['fft_size']
        self.denom = options['denom']
        self.multi = options['multi']
        self.radius = options['radius']
        self.holonomizer = Holonomizer(self.mfld_name, order=self.fft_size,
                                       radius=self.radius)
        if self.gluing_form:
            vals = [array(x) for x in self.holonomizer.R_longitude_holos]
        else:
            vals = [array(x) for x in self.holonomizer.R_longitude_evs]
        if multi == False:
            self.multiplicities, vals = self.demultiply(vals)
        self.sampled_coeffs = self.symmetric_funcs(vals)
        if self.denom:
            M = array(self.holonomizer.R_circle)
            exec('D = %s'%self.denom)
            self.raw_coeffs = array([ifft(x*D) for x in self.sampled_coeffs])
        else:
            self.raw_coeffs = array([ifft(x) for x in self.sampled_coeffs])
        self.shift = self.find_shift(self.raw_coeffs)
        if self.shift is None:
            print 'Coefficients seem to be wrapping.  A larger fft size might help.'
            return
        renorm = self.radius**(-array(range(self.fft_size - self.shift)
                                      + range(-self.shift, 0)))
        self.float_coeffs = renorm*self.raw_coeffs
        self.int_coeffs = array([map(round, x.real) for x in self.float_coeffs])
        self.height = max([max(abs(x)) for x in self.int_coeffs])
        if self.height > float(2**52):
            print "Coefficients overflowed."
        C = self.int_coeffs.transpose()
        coefficient_array =  take(C, arange(len(C))-self.shift, axis=0)
        rows, cols = coefficient_array.shape
        while rows:
            if max(abs(coefficient_array[rows-1])) > 0:
                break
            rows -= 1
        self.coefficients = coefficient_array[:rows]
        self.newton_polygon = NewtonPolygon(self.coefficients, self.gluing_form)
        self.noise = [max(abs(self.float_coeffs[i] - self.int_coeffs[i])) for
                      i in range(len(self.float_coeffs))]
        print "Noise levels: "
        for level in self.noise:
            print level
            
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
            sdr = []
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

    def find_shift(self, raw_coeffs):
       N = self.fft_size
       if N%2 == 0:
           renorm = self.radius**(-array(range(1+N/2)+range(1-N/2, 0)))
       else:
           renorm = self.radius**(-array(range(1+N/2)+range(-(N/2), 0)))
       coeffs = raw_coeffs*renorm
       maxes = [max(abs(row.real)) for row in coeffs.transpose()]
       for n in range(N):
           if maxes[-n] < 0.01:
               return n-1
    
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
            hintfile.write('"radius" : %f,\n'%self.lift.radius)
            hintfile.write('"fft_size" : %d,\n'%self.lift.fft_size)
            if not self.lift.multi:
                hintfile.write('"multi" : %s,\n'%self.lift.multi)
            if self.denom:
                hintfile.write('"denom" : "%s",\n'%self.denom)
            hintfile.write('}\n')
            hintfile.close()
            
    def boundary_slopes(self):
        print self.newton_polygon.slopes
        
    def show_lifts(self):
        # broken
        self.lift.plot()

    def show_newton(self, text=False):
        V = Polyview(self.coefficients)
        V.show_sides()
        if text:
            V.show_text()
        else:
            V.show_dots()

    def show_volumes(self):
        Plot(self.lift.volumes)

    def show_outer_volumes(self):
        Plot(self.lift.outer_volumes)

    def show_logabsL(self):
        Plot([log(abs(x)) for x in self.lift.Lsamples])

    def verify(self):
        result = True
        sign = None
        print 'Checking noise level ...',
        print max(self.noise)
        if max(self.noise) > .3:
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
    
    def tighten(self, T=1.0):
        self.holonomizer.tighten(T)
        roots = [array(x) for x in self.holonomizer.T_longitude_evs]
        self.T_sampled_coeffs = self.symmetric_funcs(roots)
        self.T_raw_coeffs = array([ifft(x) for x in self.T_sampled_coeffs])

class Slope:
      def __init__(self, xy, gluing_form=True):
            x, y = xy
            if gluing_form == False:
                  y *= 2
            if x == 0:
                  if y == 0:
                        raise ValueError, "gcd(0,0) is undefined."
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
      def __init__(self, coeff_array, gluing_form):
          self.rows, self.columns = coeff_array.shape
          self.coefficients = coeff_array
          self.gluing_form = gluing_form
          self.slopes=[]
          self.find_vertices()
            
      def find_vertices(self):
            tops = []
            bottoms = []
            for j in range(self.columns):
                  nonzero = where(self.coefficients[:,j],1,0).tolist()
                  try:
                        bottom = nonzero.index(1)
                        nonzero.reverse()
                        top = len(nonzero) - nonzero.index(1) - 1
                  except ValueError:
                        print 'Failed to construct Newton Polygon'
                        return
                  tops.append(top)
                  bottoms.append(bottom)
            if tops[0] != bottoms[0]:
                  self.slopes.append(Slope((0,1)))
            self.top_vertices = [(tops[0],0)]
            last_vertex = max = 0
            while last_vertex < self.columns - 1:
                  slope = (0,-1)
                  x,y = self.top_vertices[-1]
                  for j in range(last_vertex+1, self.columns):
                        if tops[j] == None:
                            continue
                        # Why does this sometimes throw an exception?
                        newslope =  (j - y, tops[j] - x)
                        if newslope[1]*slope[0] >= newslope[0]*slope[1]:
                              max = j
                              slope=newslope
                  self.top_vertices.append((tops[max], max))
                  s = Slope(slope, self.gluing_form)
                  if not s in self.slopes:
                        self.slopes.append(s)
                  last_vertex = max
            self.bottom_vertices = [(bottoms[0],0)]
            last_vertex = min = 0
            while last_vertex < self.columns - 1:
                  slope = (0,1)
                  x,y = self.bottom_vertices[-1]
                  for j in range(last_vertex+1, self.columns):
                        if bottoms[j] == None:
                              continue
                        newslope =  (j - y, bottoms[j] - x)
                        if newslope[1]*slope[0] <= newslope[0]*slope[1]:
                              min = j
                              slope=newslope
                  self.bottom_vertices.append((bottoms[min], min))
                  s = Slope(slope, self.gluing_form)
                  if not s in self.slopes:
                        self.slopes.append(s)
                  last_vertex = min
            self.bottom_vertices.reverse()

class Plot:
    """
    Uses gnuplot to plot a vector or list of vectors.
    Assumes that all vectors in the list are the same type (Float or Complex)
    Prompts for which ones to show.
    """
    def __init__(self, data, quiet=True, commands=''):
        self.quiet = quiet
        self.commands = commands
        if isinstance(data[0], list) or isinstance(data[0], ndarray):
            self.data = data
        else:
            self.data = [data]
        self.type = type(self.data[0][0])
        self.gnuplot = Popen(['gnuplot', '-geometry 800x720+200+0'],
                             shell=True,
                             stdin=PIPE)
        if len(self.data) > 1:
            self.show_plots()
        else:
            self.create_plot([0])
            time.sleep(1)
        self.gnuplot.terminate()
        
    def __repr__(self):
        return ''
    
    def create_plot(self, funcs):
        spec = []
        for n in funcs:
            spec.append('"-" t "%d" w lines'%n)
        gnuplot_input = self.commands + 'plot ' + ', '.join(spec) + '\n'
        if self.type == complex128 or self.type == big_complex or self.type == complex:
            for n in funcs:
                gnuplot_input += '\n'.join([
                    '%f %f'%(point.real, point.imag) if point is not None else ''
                    for point in self.data[n]] + ['e\n']) 
        elif self.type == float64 or self.type == big_float or self.type == float:
            for n in funcs:
                gnuplot_input += '\n'.join(
                    ['%f'%point for point in self.data[n]] + ['e\n']) 
        else:
            print self.type
            print self.data[0]
            raise ValueError, "Data must consist of vectors of real or complex numbers."
        self.gnuplot.stdin.write(gnuplot_input)
        
    def show_plots(self):
        if not self.quiet:
            print 'There are %d functions.'%len(self.data)
            print 'Which ones do you want to see?'
        else:
            self.create_plot( range(len(self.data)) )
        while 1:
            try:
                stuff = raw_input('plot> ')
                items = stuff.split()
                if len(items) and items[0] == 'all':
                    funcs = range(len(self.data))
                else:
                    funcs = [int(item)%len(self.data) for item in items]
                if len(funcs) == 0:
                    break
            except ValueError:
                break
            print funcs
            self.create_plot(funcs)
        return

class Polyview(NewtonPolygon):
      def __init__(self, coeff_array, gluing_form=True, scale=None, margin=50):
          self.rows, self.columns = coeff_array.shape
          self.gluing_form = gluing_form
          if scale == None:
                scale = 600/max(coeff_array.shape)
          self.scale = scale
          self.margin = margin
          self.width = (self.columns - 1)*self.scale + 2*self.margin
          self.height = (self.rows - 1)*self.scale + 2*self.margin
          self.window = Tkinter.Tk()
          self.window.wm_geometry('+400+20')
          self.coefficients = coeff_array
          self.slopes=[]
          self.find_vertices()
          self.canvas = Tkinter.Canvas(self.window,
                                  bg='white',
                                  height=self.height,
                                  width=self.width)
          self.canvas.pack()
          self.font = ('Helvetica','16','bold')
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
          return (self.margin+j*self.scale,
                  self.height - self.margin - i*self.scale)
      
      def show_dots(self):
          r = 1 + self.scale/20
          for i in range(self.rows):
              for j in range(self.columns):
                  if self.coefficients[i][j] != 0:
                      x,y = self.point((i,j))
                      self.dots.append(self.canvas.create_oval(
                          x-r, y-r, x+r, y+r, fill='black'))

      def erase_dots(self):
          for dot in self.dots:
              self.canvas.delete(dot)
          self.dots = []

      def show_text(self):
          for i in range(self.rows):
              for j in range(self.columns):
                  if self.coefficients[i][j] == 0:
                      continue
                  x,y = self.point((i,j))
                  self.text.append(self.canvas.create_text(
                          x,y,
                          text=str(self.coefficients[i][j]),
                          font=self.font,
                          anchor='c'))
      def erase_text(self):
            for coeff in self.text:
                  self.canvas.delete(coeff)
            self.text=[]

      def show_sides(self):
            r = 2 + self.scale/20
            first = self.top_vertices[0]
            x1, y1 = self.point(first)
            vertices = self.top_vertices + self.bottom_vertices + [first]
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

winding = lambda x : (sum(log(x[1:]/x[:-1]).imag) + log(x[0]/x[-1]).imag)/(-2*pi)

#M = Manifold('4_1')
#F = Fiber((-0.991020658402+0.133708842719j),
#          [
#           Point(array([
#            6.18394729421744E-01+5.14863122901458E-02j,
#            6.18394729421744E-01-5.14863122901458E-02j], dtype=DTYPE)),
#           Point(array([
#            -1.57365927858202E+00+3.47238981119960E-01j,
#            -1.57365927858202E+00-3.47238981119960E-01j], dtype=DTYPE))
#          ])
#begin = time.time()
#B = Holonomizer(M, F)
#print time.time()-begin
#Z = array(M.tetrahedra_shapes('rect'))
#print B.run_newton(Z, 1j)
