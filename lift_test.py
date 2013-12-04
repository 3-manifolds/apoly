import IPython.lib.inputhook as ih
ih.clear_inputhook()

from apoly import *
from lift_picture import *
#V = PECharVariety('m071', radius=1.04)
V = PECharVariety('m016', radius=1.02)
L = SL2RLifter(V)
P = L.show()
