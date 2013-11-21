from polish_reps import PSL2CRepOf3ManifoldGroup, PSL2RRepOf3ManifoldGroup
import snappy

sample_data = [
               ('m004(3,2)', [0.48886560625734599, 0.25766090533555303]), 
               ('m030(1,3)', [-3.0413531153050002, -0.39791078825716603, -1.0980127709735299, -1.5128763968640999]),
               ('m032(5,2)', [1.4397193529227299, -0.37244727986758702, 0.67639598574213, 0.249420323040287]),
               ('m032(7,1)', [-0.53149787942290505, 1.88940106877921, -0.367853983707157, 0.87494072873215301]),
               ('m034(-2,3)', [0.64486997829249104, 5.2944764441563397, 0.255973524222671, -0.044117577967982401]),
               ('m071(-4, 1)', [(0.55948912617811042+0.23400162016416198j), (0.67650039034242671-1.1474257693865324e-17j), (1.7704952148495239-0.94049607714155237j), (2.643772147889158+1.5020248611456067e-15j), (0.2114871099765474-2.367202429557164e-16j)])]


def sample_rep(index, precision):
    name, shapes = sample_data[index]
    M = snappy.Manifold(name)
    return PSL2CRepOf3ManifoldGroup(M, shapes, precision)
    

def basic_test():
    for i in range(5):
        rho = sample_rep(i, 10000)
        rho_real = PSL2RRepOf3ManifoldGroup(rho)
        print rho.manifold, rho_real.representation_lifts()

#basic_test()

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
    Given loxodromic matrices A and B return a basis of C^2 consisting of
    one eigenvector from each. 
    """
    basis = [ (a, b) for a in eigenvectors(A) for b in eigenvectors(B) ]
    return matrix(min(basis, key=lambda (a,b) : abs(a*b))).transpose()

def conjugate_to_PSL2R(A, B):
    C = eigenbasis(A, B)
    BB = (C**-1*A*C)
    t = 1/BB[1,1]
    C = C * matrix(B.base_ring(), [[1, 0], [0, t]])
    return (C**-1*A*C).change_ring(CDF), (C**-1*B*C).change_ring(CDF)











   



  
