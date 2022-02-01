import numpy as np
import ldpc.protograph as pt
from bposd.css import css_code
from bposd.stab import stab_code

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

    @property
    def protograph(self):
        px=pt.vstack([pt.zeros(self.hz_proto.shape),self.hx_proto])
        pz=pt.vstack([self.hz_proto,pt.zeros(self.hx_proto.shape)])
        return pt.hstack([px,pz])


class bias_tailored_lifted_product(stab_code):

    def __init__(self,lift_parameter,a,b=None):

        lhgp=lifted_hgp(lift_parameter,a,b)
        
        #Hadamard rotation
        temp1=pt.hstack([pt.zeros(lhgp.hx1_proto.shape),lhgp.hz2_proto])
        temp2=pt.hstack([lhgp.hx1_proto,pt.zeros(lhgp.hz2_proto.shape)])
        self.hx_proto=pt.vstack([temp1,temp2])
        temp1=pt.hstack([lhgp.hz1_proto,pt.zeros(lhgp.hx2_proto.shape)])
        temp2=pt.hstack([pt.zeros(lhgp.hz1_proto.shape),lhgp.hx2_proto])
        self.hz_proto=pt.vstack([temp1,temp2])

        super().__init__(self.hx_proto.to_binary(lift_parameter),self.hz_proto.to_binary(lift_parameter))

    @property
    def protograph(self):
        px=pt.vstack([pt.zeros(self.hz_proto.shape),self.hx_proto])
        pz=pt.vstack([self.hz_proto,pt.zeros(self.hx_proto.shape)])
        return pt.hstack([px,pz])





    def bias_tailor(self):
        hx_top=pt.hstack([pt.zeros(self.hx1_proto.shape),self.hz2_proto])
        hx_bottom=pt.hstack([self.hx1_proto,pt.zeros(self.hz2_proto.shape)])
        hx=self.vstack([hx_top,hx_bottom])


