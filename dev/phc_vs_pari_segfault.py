# After importing phc, a pari exception becomes a segfault

from sage.all import pari
import phc
pari([1,2])*pari([1,2])
