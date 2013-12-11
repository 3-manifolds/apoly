import euler, real_reps
import cPickle as pickle
import bz2
from apoly import *
from lift_picture import *
from lspace_slopes import compute_L_space_range as Lcone
from data import cusped_manifold_dict
import snappy
from snappy.snap.nsagetools import MapToFreeAbelianization, homological_longitude
from matplotlib.patches import Polygon

def draw_line(L, curve_on_torus, **kwargs):
    fig = L.plot.figure
    ax = fig.axis
    x = ax.get_xlim()[1]
    y = ax.get_ylim()[1]
    a, b = curve_on_torus
    if b != 0:
        ax.plot( (0, x), (0, -a*x/b), **kwargs)
    else:
        ax.plot( (0, 0), (0, -a*5), **kwargs)
    #if a != 0:
    #    ax.plot( (0, -b*y/a), (0, y), **kwargs )
    fig.draw()

def XXcone(L, C):
    u, v = C.gens()
    draw_line(L, 10000*u + v, color='black')
    draw_line(L, 10000*v + u, color='black')
    for i in range(1, 5):
        for j in range(1, 5):
            draw_line(L, i*u + j*v, color='red')
    L.plot.figure.draw()

def cone(L, C):
    #XXcone(L,C)
    fig = L.plot.figure
    ax = fig.axis
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    points = [array([0.0, 0.0])]
    A = array(C.gens())
    det = A[0,0]*A[1,1] - A[0,1]*A[1,0]
    for n, curve in enumerate(A):
        a, b = curve
        if b != 0:
            points.append(array([xmax, -a*xmax/b]))
        else:
            if a > 0:
                y = ymax if n == 0 else ymin
            else:
                y = ymin if n == 0 else ymax
#            y = 100*a if n == 0 else -100*a
            points.append(array([0, y]))
    points.insert(2, array([xmax, y]))
    vertices = array(points)
    p = Polygon(vertices, color='red', alpha=0.2, fill=True)
    L.plot.figure.axis.add_artist(p)
    L.plot.figure.draw()

def make_lifter(name):
    V = PECharVariety(name, radius=1.04)
    L = SL2RLifter(V)
    return L

def save_data(name):
    L = make_lifter(name)
    L.SL2R_rep_arcs=[]
    L.holonomizer=None
    L.SL2R_arcs = [ [ (p, list(v)) for p, v in arc] for arc in  L.SL2R_arcs]
    file = bz2.BZ2File('lifters/' + name + '.obj.bz2', 'w')
    pickle.dump(L, file)

def load_data(name):
    """
    Returns V, L.
    """
    return pickle.load(bz2.BZ2File('lifters/' + name + '.obj.bz2'))
    
def plot_data(name):
    V = PECharVariety(name, radius=1.04)
    V.show()
    L = SL2RLifter(V)
    P = L.show()
    C = Lcone(name)
    cone(L, Lcone(name))
    draw_line(L, homological_longitude(snappy.Manifold(name)), color='green')
    return P, L


def add_shifts(L, num_shifts):
    M = snappy.Manifold(L.manifold_name)
    G = M.fundamental_group()
    ab = MapToFreeAbelianization(G)
    m, l = G.peripheral_curves()[0]
    a, b = ab(m)[0], ab(l)[0]
    if a < 0:
        a, b = -a, -b
    init_arcs = L.translation_arcs
    final_arcs = init_arcs[:]
    for i in range(1, num_shifts + 1):
        for arc in init_arcs:
            final_arcs.append( [ (x+i*a, y+i*b) for x, y in arc ])
    L.translation_arcs = final_arcs

def in_homological_coordinates(L, num_shifts):
    M = snappy.Manifold(L.manifold_name)
    l = homological_longitude(M)
    assert abs(l[1]) == 1
    C = matrix(ZZ, [[1,0], l]).transpose()
    
def quick_draw(name):
    L = load_data(name)
    #add_shifts(L, 2)
    L.show()
    cone(L, Lcone(name))
    longitude = cusped_manifold_dict[name].longitude
    draw_line(L, longitude, color='green')
    return L

L = load_data('m016')
