from ldpc.codes import ring_code, rep_code
import ldpc.protograph as pt
from bposd import hgp
from lifted_hgp import lifted_hgp
import numpy as np

proto_a=pt.array([
        [(36), (), (), (), (), (0), (9)],
        [(9), (36), (), (), (), (), (0)],
        [(0), (9), (36), (), (), (), ()],
        [(), (0), (9), (36), (), (), ()],
        [(), (), (0), (9), (36), (), ()],
        [(), (), (), (0), (9), (36), ()],
        [(), (), (), (), (0), (9), (36)]
    ])
proto_b=pt.array([[(0, 1, 6)]])

qcode=lifted_hgp(lift_parameter=63,a=proto_b,b=proto_a)
qcode.test()

code_name="lifted_product_[[882,24,24]]"

for element in ["hx","hz","lx","lz"]:
    np.savetxt(f"parity_check_matrices/{code_name}_{element}.txt",eval(f"qcode.{element}"))

proto_a=pt.array([
        [(0), (11), (7), (12)],
        [(1), (8), (1), (8)],
        [(11), (0), (4), (8)],
        [(6), (2), (4), (12)]])

qcode=lifted_hgp(lift_parameter=13,a=proto_a)

qcode.test()

code_name="lifted_product_[[416,18,20]]"

for element in ["hx","hz","lx","lz"]:
    np.savetxt(f"parity_check_matrices/{code_name}_{element}.txt",eval(f"qcode.{element}"))