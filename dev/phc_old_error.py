import phc

error_message = """
Exception raised in Solve...
For the matrix A : 
 1 0 1 3 2 0 0 0 0 0 0 0
 0 0 2 2 2 0 0 1 0 0 0 0
 1 0 1 3 2 10 0 0 1 0 0 0
 1 1 0 1 -2 7 0 0 0 0 0 0
 1 0 0 0 0 9 0 0 0 0 1 0
 0 2 0 0 -1 0 0 0 0 0 0 1
 0 0 0 -1 -1 14 1 0 0 0 0 0
 0 0 -1 0 0 -2 0 0 0 0 0 0
 0 0 0 -1 -1 -10 0 0 -1 0 0 0
 0 0 0 0 2 0 0 0 0 1 0 0
 0 0 0 2 0 -4 0 0 0 0 -1 0
 0 0 0 0 -1 0 0 0 0 0 0 -1
exception happens here ...

raised ADA.NUMERICS.ARGUMENT_ERROR : a-numaux.adb:87
"""

def problem():
    equations = ['1*X0^1*X2^1*X3^1*X4^1 - 1', '1*X3^1*X5^2 - 1',
             '1*Y1^1 - 1*X0^1*X1^2*X2^1',
             '1*X0^3*X1^2*X2^3*X3^1*Y4^2 - 1*Y0^1*Y2^1',
             '1*X3^2*X5^1*Y0^1*Y2^1*Y5^1 + 1*X0^2*X1^2*X2^2*Y3^2',
             '1*X2^10*X3^7*X4^9*Y0^14 - 1*Y1^2*Y2^10*Y4^4',
             'X0 + Y0 - 1', 'X1 + Y1 - 1',
             'X2 + Y2 - 1', 'X3 + Y3 - 1', 'X4 + Y4 - 1', 'X5 + Y5 - 1']

    R = phc.PolyRing(['X0', 'X1', 'X2', 'X3', 'X4', 'X5', 'Y0', 'Y1', 'Y2', 'Y3', 'Y4', 'Y5'])
    polys = [phc.PHCPoly(R, eqn) for eqn in equations]
    system = phc.PHCSystem(R, polys)
    system.solution_list()

for i in range(10):
    print i
    problem()
             
