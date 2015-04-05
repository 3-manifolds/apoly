import apoly
import snappy

M = snappy.Manifold('m016')
P = apoly.Apoly(M)
Q = apoly.Apoly(M, gluing_form=True)

