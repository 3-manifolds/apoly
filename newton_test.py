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

def track_test(V, point=0, debug=True):
    start_time = time.time()
    H = V.holonomizer
    G = GluingSystem(V.manifold)
    steps = []
    for n in range(1,128):
        if debug: print '>>>>', n
        z = H.R_fibers[n-1].shapes[point]()
        cond = G.condition(z)
        if debug: print '1/conditions: %f %f'%(1/cond[0], 1/cond[1])
        t = H.R_circle[n]
        try:
            Zn = G.track(z,t, debug)
            if debug: print norm(Zn - H.R_fibers[n].shapes[point]())
            steps.append(Zn)
        except ValueError,e:
            print e
            break
    print time.time() - start_time
    return steps

def tighten_test(V, point=0, debug=True):
    start_time = time.time()
    H = V.holonomizer
    G = GluingSystem(V.manifold)
    steps = []
    for n in range(128):
        if debug: print '>>>>', n
        z = H.R_fibers[n].shapes[point]()
        cond = G.condition(z)
        if debug: print 'condition: %f %f'%cond
        t = H.T_circle[n]
        try:
            steps.append(G.track(z, t))
        except ValueError,e:
            print e
            break
    print time.time() - start_time
    return steps

def big_track(V):
    N = len(V.holonomizer.T_fibers[0].shapes)
    start = time.time()
    for n in range(N):
        print n
        P = track_test(V, n, debug=False)
    print time.time() - start
