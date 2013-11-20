from polish_reps import PSL2CRepOf3ManifoldGroup, PSL2RRepOf3ManifoldGroup
import snappy

sample_data = [('m030(1,3)', [-3.0413531153050002, -0.39791078825716603, -1.0980127709735299, -1.5128763968640999]),
               ('m032(5,2)', [1.4397193529227299, -0.37244727986758702, 0.67639598574213, 0.249420323040287]),
               ('m032(7,1)', [-0.53149787942290505, 1.88940106877921, -0.367853983707157, 0.87494072873215301]),
               ('m034(-2,3)', [0.64486997829249104, 5.2944764441563397, 0.255973524222671, -0.044117577967982401]),
               ('m071(-4, 1)', [0.559453254962+0.233744240384j, 0.676554283091, 1.77127144543-0.939796975438j, 2.64352893406, 0.211509177188])]


def sample_rep(index, precision):
    name, shapes = sample_data[index]
    M = snappy.Manifold(name)
    return PSL2CRepOf3ManifoldGroup(M, shapes, precision)
    

def basic_test():
    for i in range(4):
        rho = sample_rep(i, 500)
        print rho, rho.is_PSL2R_rep()
        rho_real = PSL2RRepOf3ManifoldGroup(rho)
        print rho_real.representation_lifts()

basic_test()
