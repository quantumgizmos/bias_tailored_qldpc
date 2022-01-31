# Bias-tailored quantum LDPC codes

This repository cotains the code for the decoding simulations of bias-tailored quantum LDPC codes as described in arxiv:2202:xxxx. Also included are various examples showing how our code base can be used to construct lifted product codes.

## Setup
The code in this repository requires functions from [`LDPC`](https://github.com/quantumgizmos/ldpc) and [`BPOSD`](https://github.com/quantumgizmos/bp_osd). To install, run the following command

```
pip install -U ldpc bposd
```

# Classical error correction 
Section 2.1 in arXiv:2202:xxxx



A binary linear code is defined by its parity check matrix. Eg. the parity check matrix for the four-bit (closed loop) repetition code is given by


```python
import numpy as np
from ldpc.codes import ring_code

H=ring_code(4)
print(H)
```

    [[1 1 0 0]
     [0 1 1 0]
     [0 0 1 1]
     [1 0 0 1]]


The code parameters of this code are [n=4,k=1,d=4] where n is the block length, k is the number of encoded bits and d is the code distance. These parameters can be checked as follows:


```python
from ldpc.code_util import compute_code_distance
import ldpc.mod2 as mod2

n=H.shape[1] #the code block length
k=n-mod2.rank(H) #This follows from the rank-nullity theorem
d=compute_code_distance(H) #This function exhaustively computes the code distance by checking the weight of all logical operators (Not scalable!).
print(f"[{n},{k},{d}]")
```

    [4,1,4]


## Quasi-cyclic codes
Section 2.3 in arXiv:2202:xxxx

Quasicyclic codes are a type of classical linear code that are highly performant under BP decoding. The parity check matrix of a quasi-cyclic code is a block matrix where each block corresponds to a sum over permutation matrices.

### Permutation matrices

A permuation matrix is defined as the matrix obtained by shifting an `LxL` matrix by `t` positions to the right. After `t=L` shifts, the identity is recovered. eg.


```python
from ldpc import protograph as pt
for i in range(4):
    print(f"t={i}")
    print(pt.permutation_matrix(3,i))
```

    t=0
    [[1 0 0]
     [0 1 0]
     [0 0 1]]
    t=1
    [[0 1 0]
     [0 0 1]
     [1 0 0]]
    t=2
    [[0 0 1]
     [1 0 0]
     [0 1 0]]
    t=3
    [[1 0 0]
     [0 1 0]
     [0 0 1]]


### The ring of circulants

The set of `L` distinct permutation matrices forms a basis of a vector space referred to as the ring of circulants. This group is isomorphic to the polynomial equation `x^L-1=0`, providing a compact algebra with which to represent it. We can define a ring element as follows:


```python
ring_element=pt.ring_of_circulants_f2((0,1))
print(ring_element)
```

    (0,1)


The above ring element can be mapped to binary representation (of size L) as follows:


```python
L=3
ring_element.to_binary(3)
```




    array([[1, 1, 0],
           [0, 1, 1],
           [1, 0, 1]])



We see that the ring element corresponds to sum of permutation matrices with shift parameters `t=0` and `t=1`. In general, all computations on permutation matrices involving addition and multiplication mod2 can be represented using the ring algebra. For example:

#### Matrix version


```python
L=5

a=pt.permutation_matrix(L,2)+pt.permutation_matrix(L,4)
b=pt.permutation_matrix(L,3)

print((a@b %2+(a+b)%2)%2)
```

    [[1 0 0 1 1]
     [1 1 0 0 1]
     [1 1 1 0 0]
     [0 1 1 1 0]
     [0 0 1 1 1]]


#### Ring version


```python
L=5

a=pt.ring_of_circulants_f2((2,4))
b=pt.ring_of_circulants_f2((3))
c= a*b + (a+b)

print(c)

print(c.to_binary(L))
```

    (5,7,2,4,3)
    [[1 0 0 1 1]
     [1 1 0 0 1]
     [1 1 1 0 0]
     [0 1 1 1 0]
     [0 0 1 1 1]]


### Protographs

A protograph is an `mxn` matrix where each element is in the ring of circulants. eg.


```python
A=pt.array([[(1,2),(0),()],
             [(0),(0,1),(1)]])

print(a)
```

    (2,4)


We can convert this to a binary parity check matrix as follows:


```python
H=A.to_binary(lift_parameter=3)
print(H)
```

    [[0 1 1 1 0 0 0 0 0]
     [1 0 1 0 1 0 0 0 0]
     [1 1 0 0 0 1 0 0 0]
     [1 0 0 1 1 0 0 1 0]
     [0 1 0 0 1 1 0 0 1]
     [0 0 1 1 0 1 1 0 0]]


The parity check matrix obtained from the binary representation of the protograph defines a quasi-cyclic code.

### The repetition code as a quasi-cyclic code

The family of length-n [n,1,n] closed-loop repetition codes can be represented as a quasi-cyclic code family with protographs of the form:


```python
A=pt.array([[(0,1)]])
H=A.to_binary(lift_parameter=4)
print(H)
```

    [[1 1 0 0]
     [0 1 1 0]
     [0 0 1 1]
     [1 0 0 1]]


# Quantum error correction

## Calderbank, Shor & Steane (CSS codes)

Section 3.2 in arXiv:2202.xxxx

CSS codes are quantum stabiliser codes with non-overlapping X- and Z- stabilisers. Using our code, you can construct a CSS code as follows:


```python
from ldpc.codes import hamming_code
H=hamming_code(3)
print(H)
```

    [[0 0 0 1 1 1 1]
     [0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1]]



```python
from bposd.css import css_code

qcode=css_code(hx=H,hz=H)

print(qcode.hx)
print(qcode.hz)
```

    [[0 0 0 1 1 1 1]
     [0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1]]
    [[0 0 0 1 1 1 1]
     [0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1]]


In the above `hx` and `hz` are parity check matrices that describe the X- and Z-stabiliser measurments respectively. The logical operators of a CSS code are obtained from our `css_code` object as follows:


```python
lx=qcode.lx #x logical operators
lz=qcode.lz #z logical operators

print(lx)
print(lz)
```

    [[1 1 1 0 0 0 0]]
    [[1 1 1 0 0 0 0]]


All stabiliser codes must have `hx` and `hz` matrices that commute. We can test that this requirement is satisfied using the `css_code.test` function. eg.


```python
qcode.test()
```

    <Unnamed CSS code>, (3,4)-[[7,1,nan]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (3,4)-[[7,1,nan]]





    True



From the above we see that the code we have defined is a valid CSS code. However, this will not be the case for arbitrary combinations of `hx` and `hz` matrices. For, example, if we define `hx` and `hz` to be repetition codes, we get the following invalid CSS code:


```python
h=ring_code(4)
qcode=css_code(hx=h,hz=h)
qcode.test()

```

    <Unnamed CSS code>, (2,2)-[[4,-2,nan]]
     -Block dimensions incorrect
     -PCMs commute hz@hx.T==0: Fail
     -PCMs commute hx@hz.T==0: Fail
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anitcommute: Fail





    False



The `qcode.test` function executed above tells us that the `hx` and `hz` matrices do not commute.

## Hypergraph product codes
Section 3.3 in arXiv:2202:xxxx

We have now seen that CSS stabiliser codes cannot take any arbitrary combination of `hx` and `hz` matrices as inputs due to the requirement that stabilisers must commute. So how do we create quantum codes from the starting point of classical codes? One solution is to use a code construction method called the hypegraph product which was first proposed by Tillich and Zemor. This defines the `hx` and `hz` components of the CSS code as follows:

```
hx=[ h1 ⊗ I, h2^T ⊗ I ]
hz=[ I ⊗ h2, h1^t ⊗ I]
```

### The toric code form the hypergraph product

The toric code can be constructed from the hypergraph product of two closed loop repetition codes.

where `h1` and `h2` are two classical *seed* codes. It is straightforward to verify that the above two parity check matrices will commute for any combination of the seed codes. The hypergraph product therefore allows us to construct a quantum code from any arbitrary pair of classical binary codes.


```python
from bposd.hgp import hgp
h1=ring_code(2)
h2=ring_code(3)

qcode=hgp(h1,h2,compute_distance=True)
qcode.test()

```

    <Unnamed CSS code>, (2,4)-[[12,2,2]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (2,4)-[[12,2,2]]





    True



### A quantum LDPC code from the hypergraph product

The hypergraph product can also be used to create quantum LDPC codes from the starting point of classical LDPC codes. For example, consider the following classical LDPC code of the length 16


```python
H=np.array([[1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0],
            [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1],
            [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
            [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]])

from ldpc.code_util import get_code_parameters

n,k,d,_,_=get_code_parameters(H)

print(f"Code parameters: [{n},{k},{d}]")
```

    Code parameters: [16,4,6]


The hypergraph product of two pairs of the above code gives the following quantum LDPC code


```python
qcode=hgp(H,H,compute_distance=True)
qcode.test()
```

    <Unnamed CSS code>, (4,7)-[[400,16,6]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (4,7)-[[400,16,6]]





    True




```python

```
