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
            exec('D = %s'%self.denom) #denom is an expression for a poly in H_M
            self.raw_coeffs = array([ifft(x*D) for x in self.sampled_coeffs])
        else:
            self.raw_coeffs = array([ifft(x) for x in self.sampled_coeffs])
        #print self.raw_coeffs.shape
        self.shift = self.find_shift()
        if self.shift is None:
            print 'Coefficients seem to be wrapping.  '
                  'A larger fft size might help.'
            return
        print 'Shift is %s.'%self.shift
        renorm = array([self.radius])**(-array(range(self.fft_size - self.shift)
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
        # If the FFT size is large enough and we have a good denominator,
        # we should not go past the middle
        realrows = rows/2 - 1
        while realrows:
            if max(abs(coefficient_array[rows-1])) > 0:
                break
            if realrows:
                realrows -= 1
        self.coefficients = coefficient_array[:realrows]
        self.noise = [max(abs(self.float_coeffs[i] - self.int_coeffs[i])) for
                      i in range(len(self.float_coeffs))]
        print "Noise levels: "
        for level in self.noise:
            print level
        if max(self.noise) > 0.1:
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

    def Xfind_shift(self, raw_coeffs):
       rows, cols = raw_coeffs.shape
       N = self.fft_size
       R = array([self.radius])
       shifts = [0]
       if N%2 == 0:
           renorm = R**(-array(range(1+N/2)+range(1-N/2, 0)))
       else:
           renorm = R**(-array(range(1+N/2)+range(-(N/2), 0)))
       coeffs = raw_coeffs*renorm
       for i in range(rows):
          for j in range(1, 1+ cols/2):
             if abs(abs(coeffs[i][-j]) - 1.) < .001:
                 shifts.append(j)
       print 'shifts: ', shifts
       return max(shifts)

    def find_shift(self, cutoff=0.1):
       N = self.fft_size
       R = self.radius
       if N%2 == 0:
           renorm = R**(-array(range(1+N/2)+range(1-N/2, 0)))
       else:
           renorm = R**(-array(range(1+N/2)+range(-(N/2), 0)))
       coeffs = (self.raw_coeffs*renorm).transpose()
       if max(abs(coeffs[N/2])) > 0.5:
              return None
       for n in range(1+N/2,N):
           if max(abs(coeffs[n])) > cutoff:
               return N - n
       return 0
    
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
        V = PolyViewer(self.newton_polygon())
        V.show_sides()
        if text:
            V.show_text()
        else:
            V.show_dots()

class ShapeRelation(list):
    def __init__(self, manifold_name):
        self.holonomizer = H = Holonomizer(manifold_name)
        H.tighten()
        self.shapes = [
            [[f.shapes[m][n] for f in H.T_fibers] for m in range(H.degree)]
            for n in range(H.dim)]
        for shape in self.shapes:
            self.append(PolyRelation(shape, H.T_circle, radius=1.0))
