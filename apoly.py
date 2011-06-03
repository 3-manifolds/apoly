from numpy import array, matrix, dot, prod, diag, transpose, zeros, ones
from numpy import log, exp, pi, sqrt
from numpy.linalg import svd, norm
from numpy import dtype, complex128, float64
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
import time, sys


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
    equality operator.  Instantiate with a sequence of complexes.
    """
    def __init__(self, Z):
        self.Z = array(Z)
        
    def __repr__(self):
        return str(self.Z)

    def __eq__(self, other):
        # NOTE: For knot 9_30 the fibers only matched to 1.9E-8
        return norm(self.Z - other.Z) < 1.0E-10

    def __xor__(self, other):
        return norm(self.Z - other.Z)
    
class Fiber:
    def __init__(self, H_meridian, system, tolerance=1.0E-06):
        """
        A fiber for the rational function [holonomy of the meridian].
        Manages a single solved PHCSytem.
        """
        self.H_meridian = H_meridian
        self.system = system
        N = system.num_variables()/3
        self.solutions = system.solution_list(tolerance=tolerance)
        # only keep the "X" variables.
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
            for z in p.Z:
                if abs(z*(1-z)) < 1.0E-10:
                    return False
        return True
            
    def details(self):
        for n, s in enumerate(self.solutions):
            print 'solution #%s:'%n
            print s

    def residuals(self):
        for n, s in enumerate(self.solutions):
            print n, s.res 

class PHCFibrator:
    """
    A factory for fibers computed by PHC.
    """
    def __init__(self, mfld_name, radius=1.02):
        self.mfld_name = mfld_name
        self.radius = radius
        self.manifold = Manifold(self.mfld_name)
        self.num_tetrahedra = N = self.manifold.num_tetrahedra()
        variables = ( ['X%s'%n for n in range(N)] +
                      ['Y%s'%n for n in range(N)] +
                      ['Z%s'%n for n in range(N)] )
        self.ring = PolyRing(variables + ['t'])
        self.basepoint = radius*exp(2*pi*1j*random())
        self.equations = self.build_equations()
        self.equations += ['X%s*Y%s*Z%s - 1'%(n,n,n) for n in range(N)]
        self.equations += ['X%s + Y%s - 1'%(n,n) for n in range(N)] 
        self.psystem = ParametrizedSystem(
            self.ring,
            't',
            [PHCPoly(self.ring, e) for e in self.equations]
            )
        #print 'Computing the starting fiber ... ',
        #begin = time.time()
        #self.base_system = self.psystem.start(self.basepoint, tolerance=1.0E-05)
        #print 'done. (%s seconds)'%(time.time() - begin)
        #self.base_fiber = Fiber(self.basepoint,
        #                         self.base_system)

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
            raise RuntimeError, 'Starting fiber contains ideal points.'
        self.degree = len(self.base_fiber)
        # pre-initialize
        self.R_fibers = range(order)
        self.fibers = range(order)
        self.dim = manifold.num_tetrahedra()
        self.rhs = []
        eqns = manifold.gluing_equations('rect')
        self.glunomials = [Glunomial(A, B, c) for A, B, c in eqns[:-3]]
        self.rhs = [1.0]*(len(eqns) - 3)
        self.perturbation = 0.1*array([exp(2*pi*1j*random())
                                          for n in xrange(self.dim - 1)] + [0])
        self.M_holo, self.L_holo = [Glunomial(A,B,c) for A,B,c in eqns[-2:]]
        self.glunomials.append(self.M_holo)
        print 'Tracking satellite ...'
        self.track_satellite()
        print 'Tightening things up ...'
        self.tighten()
        print 'Computing longitude holonomies and eigenvalues'
        try:
            self.R_longitude_holos, self.R_longitude_evs = self.longitude_data(self.R_fibers)
            self.longitude_holos, self.longitude_evs = self.longitude_data(self.fibers)
        except:
            print 'Failed'
            
    def __call__(self, Z):
        return array([F(Z) for F in self.glunomials])

    def track_satellite(self):
        arg = log(self.base_fiber.H_meridian).imag%(2*pi)
        R = self.fibrator.radius
        Darg = 2*pi/self.order
        self.base_index = base_index = int(arg/Darg)
        self.R_circle = circle = [R*exp(n*Darg*1j) for n in range(self.order)]
        H_base = circle[base_index]
        self.R_fibers[base_index] = self.base_fiber
        for n in xrange(base_index+1, self.order):
            print n,
            sys.stdout.flush()
            self.R_fibers[n] = F = self.fibrator.transport(
                self.R_fibers[n-1], circle[n])
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
        if not self.last_R_fiber == self.R_fibers[0]:
            print 'Lifts did not close up'
            print array(self.last_R_fiber.points)
            print array(self.R_fibers[0].points)

    def tighten(self):
        Darg = 2*pi/self.order
        self.circle = circle = [exp(n*Darg*1j) for n in range(self.order)]
        for n in xrange(self.order):
            print n,
            sys.stdout.flush()
            self.fibers[n] = self.fibrator.transport(self.R_fibers[n],
                                                     circle[n])
        print
        self.last_fiber = self.fibrator.transport(self.R_fibers[-1],
                                                  self.circle[0])
        if not self.last_fiber == self.fibers[0]:
            print 'Lifts did not close up'
            print array(self.last_fiber.points)
            print array(self.fibers[0].points)

    def longitude_data(self, fiber_list):    
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
            for n, H in enumerate(L, 1):
                e = sqrt(H)
                E.append(e if abs(e - E[n-1]) < abs(e + E[n-1]) else -e )
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

    def newton_step(self, Z, target):
        U, s, V = svd(self.jacobian(Z))
        D = matrix(diag(1/s))
        rhs = transpose(target - self(Z))
        return s, array(V.H*D*U.H*rhs).T[0]

    def run_newton(self, Z, target, info):
        delta = self(Z) - target
        prior_residual = norm(delta)
        info['residuals'] = residuals = []
        increments = []
        singular = []
        result = None
        N = 0
        while True:
            N += 1
            singular_values, dZ = self.newton_step(Z, target)
            min_sing = min(singular_values)
            singular.append(min_sing)
            if min_sing < SVT:
                self.crash_log.append([residuals, increments,
                                 singular_values, Z,
                                 backtracks, N])
                raise Singularity
            nextZ = Z + dZ
            residual = norm(self(nextZ) - target)
            if residual < 0.0001 and residual > prior_residual:
                residuals.append(residual)
                return Z
            norm_dZ = norm(dZ)
            for backtracks in range(1,33):
                if residual < 1.0E-15 or norm_dZ < 1.0E-6:
                    break
                else:
                    next_residual = norm(self(Z + 0.5*dZ) - target)
                    if next_residual > residual:
                        break
                    residual = next_residual
                    dZ *= 0.5
                    norm_dZ *= 0.5
                    nextZ = Z + dZ
            residuals.append(residual)
            increment = norm(dZ)
            increments.append( increment )
            if N >= MAX_STEPS or residual < 1.0E-15:
                info['increments'] = increments
                info['singular_values'] = singular
                info['num_steps'] = N
                info['backtracks'] = backtracks
                return nextZ
            prior_residual = residual
            Z = nextZ

class Plot:
    """
    Uses gnuplot to plot a vector or list of vectors.
    Assumes that all vectors in the list are the same type (Float or Complex)
    Prompts for which ones to show.
    """
    def __init__(self, data, quiet=False):
        if isinstance(data[0], list) or isinstance(data[0], array):
            self.data = data
        else:
            self.data =[data]
        self.type = type(data[0][0])
        self.gnuplot = Popen(['gnuplot', '-geometry 800x720+200+0'],
                             shell=True,
                             stdin=PIPE)
        self.show_plots(quiet)

    input = raw_input
    
    def show_plots(self, quiet):
        if not quiet:
            print 'There are %d functions.'%len(self.data)
        print 'Which ones do you want to see?'
        while 1:
            try:
                stuff = self.input('plot> ')
                items = stuff.split()
                if len(items) and items[0] == 'all':
                    list = range(len(self.data))
                else:
                    list = [int(item)%len(self.data) for item in items]
                if len(list) == 0:
                    break
            except ValueError:
                break
            print list
            spec = []
            for n in list:
                spec.append('"-" t "%d" w lines'%n)
            gnuplot_input = 'plot ' + ', '.join(spec) + '\n'
            if self.type == complex128 or self.type == big_complex:
                for n in list:
                    gnuplot_input += '\n'.join([
                        '%f %f'%(point.real, point.imag)
                        for point in self.data[n]] + ['e\n']) 
            elif self.type == float64 or self.type == big_float:
                    gnuplot_input += '\n'.join(
                        ['%f'%point for point in self.data[n]] + ['e\n']) 
            else:
                print self.type
                print self.data[0]
                print "Data must consist of vectors of real or complex numbers."
                return
            self.gnuplot.stdin.write(gnuplot_input)
        #self.gnuplot.close()
        return

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
