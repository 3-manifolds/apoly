import cPickle as pickle
import bz2, glob, os
from apoly import *
from lift_picture import *
from lspace_slopes import compute_L_space_range as Lcone
from census_data import cusped_manifold_dict
import snappy
from snappy.snap.nsagetools import MapToFreeAbelianization
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

def save_data(name, to_file=True, save_holonomizer=False):
    L = make_lifter(name)
    L.SL2R_rep_arcs=[]
    if save_holonomizer:
        L.holonomizer.base_fiber = None
        L.holonomizer.fibrator = None
        L.holonomizer.manifold = None
        L.holonomizer.hp_manifold = None
    else:
        L.holonomizer = None
    if to_file:
        file = bz2.BZ2File('lifters/' + name + '.obj.bz2', 'w')
        pickle.dump(L, file, 2)
    else:
        return bz2.compress(pickle.dumps(L, 2))

def precomputed_lifters():
    return [match[8:-8] for match in glob.glob('lifters/*.obj.bz2')]
    
def load_data(name):
    """
    Returns V, L.
    """
    return pickle.load(bz2.BZ2File('lifters/' + name + '.obj.bz2'))
    
def plot_data(name, radius=1.04):
    V = PECharVariety(name, radius=radius)
    V.show()
    L = SL2RLifter(V)
    P = L.show()
    draw_line(L, snappy.Manifold(name).homological_longitude(), color='green')
    try:
        C = Lcone(name)
        cone(L, C)
    except:
        pass
    if cusped_manifold_dict.has_key(name):
        for slope in cusped_manifold_dict[name].L_space_fillings:
            draw_line(L, slope, color='red')
        for slope in cusped_manifold_dict[name].non_L_space_fillings:
            draw_line(L, slope, color='#A9F5F2')
    return L


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
    if isinstance(name, str):
        L = load_data(name)
    else:
        L = name
        name = name.manifold.name()
        
    L.show()
    try:
        cone(L, Lcone(name))
    except:
        pass
    longitude = cusped_manifold_dict[name].longitude
    draw_line(L, longitude, color='green')
    for slope in cusped_manifold_dict[name].L_space_fillings:
        draw_line(L, slope, color='red')
    for slope in cusped_manifold_dict[name].non_L_space_fillings:
        draw_line(L, slope, color='#A9F5F2')
    return L

def save_all_plots():    
    from matplotlib.backends.backend_pdf import PdfPages
    pdf_pages = PdfPages('all_plots.pdf')
    for name in precomputed_lifters():
        L = quick_draw(name)
        L.plot.figure.figure.savefig(pdf_pages, format='pdf')
    pdf_pages.close()

def save_all_plots_new():    
    from matplotlib.backends.backend_pdf import PdfPages
    import taskdb2

    db = taskdb2.TaskDatabase('PEChar')
    df = db.dataframe('done=1')
    pdf_pages = PdfPages('all_plots_new.pdf')
    for i, row in df.iterrows():
        L = pickle.loads(bz2.decompress(row.lifter))
        if len(L.translation_arcs):
            try:
                L = quick_draw(L)
                L.plot.figure.figure.savefig(pdf_pages, format='pdf')
            except:
                print "Issue with " + row['name']
    pdf_pages.close()
    
#save_all_plots()
    

#L = load_data('m016')


