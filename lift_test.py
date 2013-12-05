try:
    import IPython.lib.inputhook as ih
    ih.clear_inputhook()
except:
    pass
import euler, real_reps
import cPickle as pickle
import bz2
from apoly import *
from lift_picture import *
from lspace_slopes import compute_L_space_range as Lcone
from lspace_data import longitudes



def draw_line(L, curve_on_torus, **kwargs):
    ax = L.plot.figure.axis
    x = ax.get_xlim()[1]
    y = ax.get_ylim()[1]
    a, b = curve_on_torus
    if b != 0:
        ax.plot( (0, x), (0, -a*x/b), **kwargs)
    else:
        ax.plot( (0, x), (0, 0), **kwargs)
    #if a != 0:
    #    ax.plot( (0, -b*y/a), (0, y), **kwargs )
        

def cone(L, C):
    u, v = C.gens()
    draw_line(L, u, color='black')
    draw_line(L, v, color='black')
    draw_line(L, u + v, color='red')
    draw_line(L, 2*u + v, color='red')
    draw_line(L, 3*u + v, color='red')
    draw_line(L, u + 2*v, color='red')
    draw_line(L, u + 3*v, color='red')
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
    return P, L


def quick_draw(name):
    L = load_data(name)
    L.show()
    cone(L, Lcone(name))
    draw_line(L, longitudes[name], color='green')
    return L
