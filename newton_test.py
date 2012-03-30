def T_test(V):
    H = V.holonomizer
    G = GluingSystem(V.manifold)
    for n in range(1,128):
        print '>>>>', n
        z = H.T_fibers[n-1].shapes[0]()
        print 'condition: %f'%G.condition(z)
        t = H.T_circle[n]
        zz = G.newton1(z, t) 

def R_test(V):
    H = V.holonomizer
    G = GluingSystem(V.manifold)
    for n in range(1,128):
        print '>>>>', n
        z = H.R_fibers[n-1].shapes[0]()
        print 'condition: %f'%G.condition(z)
        t = H.R_circle[n]
        zz = G.newton2(z, t) 

def track_test(V):
    H = V.holonomizer
    G = GluingSystem(V.manifold)
    for n in range(1,128):
        print '>>>>', n
        z = H.R_fibers[n-1].shapes[0]()
        print 'condition: %f'%G.condition(z)
        t = H.R_circle[n]
        zz = G.track(z, t) 
