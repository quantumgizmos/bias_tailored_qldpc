import numpy as np
import protograph as pt
from bposd.css import css_code

def I(n):
    return pt.identity(n)

class lifted_hgp(css_code):

    def __init__(self,lift_parameter,a,b=None):

        '''
        Generates the lifted hypergraph product of the protographs a and b
        '''
        self.a=a

        self.a_m,self.a_n=self.a.shape

        if b is None:
            self.b=pt.copy(self.a)
        else:
            self.b=b
        
        self.b_m,self.b_n=self.b.shape

        self.hx1_proto=pt.kron(self.a,I(self.b_n))
        self.hx2_proto=pt.kron(I(self.a_m),self.b.T)
        self.hx_proto=pt.hstack([self.hx1_proto,self.hx2_proto])

        self.hz1_proto=pt.kron(I(self.a_n),self.b)
        self.hz2_proto=pt.kron(self.a.T,I(self.b_m))
        self.hz_proto=pt.hstack([self.hz1_proto,self.hz2_proto])

        super().__init__(self.hx_proto.to_binary(lift_parameter),self.hz_proto.to_binary(lift_parameter))