
# Bias-tailored quantum LDPC codes

This repository contains the code for the decoding simulations of bias-tailored quantum LDPC codes as described in arxiv:2202:xxxx. Also included are various examples showing how our code base can be used to construct lifted product codes.

- [Paste Your Document In Here](#paste-your-document-in-here)
  * [And a table of contents](#and-a-table-of-contents)
  * [On the right](#on-the-right)
- [Bias-tailored quantum LDPC codes](#bias-tailored-quantum-ldpc-codes)
  * [Setup](#setup)
- [Classical error correction](#classical-error-correction)
  * [Quasi-cyclic codes](#quasi-cyclic-codes)
    + [Permutation matrices](#permutation-matrices)
    + [The ring of circulants](#the-ring-of-circulants)
      - [Matrix version](#matrix-version)
      - [Ring version](#ring-version)
    + [Protographs](#protographs)
    + [The repetition code as a quasi-cyclic code](#the-repetition-code-as-a-quasi-cyclic-code)
- [Quantum error correction](#quantum-error-correction)
  * [Calderbank, Shor & Steane (CSS codes)](#calderbank--shor---steane--css-codes-)
  * [Hypergraph product codes](#hypergraph-product-codes)
    + [The toric code form the hypergraph product](#the-toric-code-form-the-hypergraph-product)
    + [A quantum LDPC code from the hypergraph product](#a-quantum-ldpc-code-from-the-hypergraph-product)
  * [Lifted product codes](#lifted-product-codes)
- [Bias-tailoring](#bias-tailoring)
  * [The CSS twisted toric code](#the-css-twisted-toric-code)
    + [The CSS twisted toric code under infinite-bias](#the-css-twisted-toric-code-under-infinite-bias)
  * [The XZZX twisted toric code](#the-xzzx-twisted-toric-code)
    + [The XZZX twisted toric code under infinite bias](#the-xzzx-twisted-toric-code-under-infinite-bias)
    + [General construction for XZZX twisted toric codes](#general-construction-for-xzzx-twisted-toric-codes)
  * [Bias-tailored LDPC codes](#bias-tailored-ldpc-codes)
- [BP+OSD decoding of bias-tailored LDPC codes](#bp-osd-decoding-of-bias-tailored-ldpc-codes)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


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



## Lifted product codes

The lifted product is an extension of the hypergraph product that allows quantum codes to be constructed from pairs of classical quasi-cyclic codes. The lifted product was first proposed by Pantaleev and Kalachev and has been proved to yield asympotically good quantum LDPC codes with linear distance-to-block length scaling. In a lifted product code, the `hx` and `hz` CSS components are constructed from two protographs defined as follows:

```
ax=[ a1 ⊗ E, a2^T ⊗ E ]
az=[ E ⊗ a2, a1^t ⊗ E]
```

where `a1` and `a2` are the protographs of classical quasi-cyclic codes and `E` is the identity protograph element. The `hx` and `hz` components are obtained by mapping `ax` and `az` to their binary representation. Similiar to the hypergraph product, the lifted product allows a quantum code to be constructed from any pair of protographs. As an example, consider the protograph below


```python
a1=pt.array([
        [(0), (11), (7), (12)],
        [(1), (8), (1), (8)],
        [(11), (0), (4), (8)],
        [(6), (2), (4), (12)]])

a1
```




    protograph.array([[(0),(11),(7),(12)],[(1),(8),(1),(8)],[(11),(0),(4),(8)],[(6),(2),(4),(12)]])



Mapping the above protograph to a binary with lift parameter L=13 gives the following classical code 


```python
H=a1.to_binary(lift_parameter=13)
n,k,d,_,_=get_code_parameters(H)
print(f"Code parameters: [{n},{k},{d}]")
```

    Code parameters: [52,3,26]


Now if we take the hypergraph product of this binary matrix, we get the following CSS code


```python
# qcode=hgp(H,H,compute_distance=True) #this will take a while
# qcode.test()
```

The rate of this code is:


```python
rate=qcode.K/qcode.N
rate
```




    0.04



It is clear that the rate of this code is quite small! In the lifted product, we instead define the quantum code as a hypergraph product of two protographs. This requires performing all tensor product operations over the ring algebra, but results in a much more compact code.


```python
from lifted_hgp import lifted_hgp
quantum_protograph_code=lifted_hgp(lift_parameter=13,a=a1,b=a1)
print(quantum_protograph_code.hz_proto)
```

    [[(0) (11) (7) (12) () () () () () () () () () () () () (0) () () () (-1) () () () (-11) () () () (-6) () () ()]
     [(1) (8) (1) (8) () () () () () () () () () () () () () (0) () () () (-1) () () () (-11) () () () (-6) () ()]
     [(11) (0) (4) (8) () () () () () () () () () () () () () () (0) () () () (-1) () () () (-11) () () () (-6) ()]
     [(6) (2) (4) (12) () () () () () () () () () () () () () () () (0) () () () (-1) () () () (-11) () () () (-6)]
     [() () () () (0) (11) (7) (12) () () () () () () () () (-11) () () () (-8) () () () (0) () () () (-2) () () ()]
     [() () () () (1) (8) (1) (8) () () () () () () () () () (-11) () () () (-8) () () () (0) () () () (-2) () ()]
     [() () () () (11) (0) (4) (8) () () () () () () () () () () (-11) () () () (-8) () () () (0) () () () (-2) ()]
     [() () () () (6) (2) (4) (12) () () () () () () () () () () () (-11) () () () (-8) () () () (0) () () () (-2)]
     [() () () () () () () () (0) (11) (7) (12) () () () () (-7) () () () (-1) () () () (-4) () () () (-4) () () ()]
     [() () () () () () () () (1) (8) (1) (8) () () () () () (-7) () () () (-1) () () () (-4) () () () (-4) () ()]
     [() () () () () () () () (11) (0) (4) (8) () () () () () () (-7) () () () (-1) () () () (-4) () () () (-4) ()]
     [() () () () () () () () (6) (2) (4) (12) () () () () () () () (-7) () () () (-1) () () () (-4) () () () (-4)]
     [() () () () () () () () () () () () (0) (11) (7) (12) (-12) () () () (-8) () () () (-8) () () () (-12) () () ()]
     [() () () () () () () () () () () () (1) (8) (1) (8) () (-12) () () () (-8) () () () (-8) () () () (-12) () ()]
     [() () () () () () () () () () () () (11) (0) (4) (8) () () (-12) () () () (-8) () () () (-8) () () () (-12) ()]
     [() () () () () () () () () () () () (6) (2) (4) (12) () () () (-12) () () () (-8) () () () (-8) () () () (-12)]]


Mapping the above protograph to binary gives us the `hx` component of the CSS code.


```python
hx=quantum_protograph_code.hx_proto.to_binary(lift_parameter=13)
hx
```




    array([[1, 0, 0, ..., 0, 0, 0],
           [0, 1, 0, ..., 0, 0, 0],
           [0, 0, 1, ..., 0, 0, 0],
           ...,
           [0, 0, 0, ..., 0, 1, 0],
           [0, 0, 0, ..., 0, 0, 1],
           [0, 0, 0, ..., 0, 0, 0]])



Similarily for `hz`


```python
hz=quantum_protograph_code.hz_proto.to_binary(lift_parameter=13)
hz
```




    array([[1, 0, 0, ..., 0, 0, 0],
           [0, 1, 0, ..., 0, 0, 0],
           [0, 0, 1, ..., 0, 0, 0],
           ...,
           [0, 0, 0, ..., 0, 1, 0],
           [0, 0, 0, ..., 0, 0, 1],
           [0, 0, 0, ..., 0, 0, 0]])



Now, constructing a CSS code from this pair `hx` and `hz` gives the following


```python
qcode=css_code(hx,hz)
qcode.test()
```

    <Unnamed CSS code>, (4,8)-[[416,18,nan]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (4,8)-[[416,18,nan]]





    True



Note, we can skip the above step and obtain the CSS code directly from the `lifted_hgp` object 


```python
qcode=lifted_hgp(lift_parameter=13,a=a1,b=a1)
qcode.test()

```

    <Unnamed CSS code>, (4,8)-[[416,18,nan]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (4,8)-[[416,18,nan]]





    True



We have now seen two examples of how a quantum code can be constructed from the [52,3,6] quasi-cyclic code:
- Taking the hypergraph product of the code's binary parity check matrix yields a CSS code with parameters [[[5408,18,26]]]
- Taking the lifted product of the code's protgraph yields a CSS code with [[416,18,d~20]]

Note that the distance of `d~20` for the lifted product is an estimate based on the lowest weight numerically observed logical operator (more on that below). Whilst the distance of the lifted product is less than that for the hypergrap product (`d~20` vs. `d=26`), the lifted product has a much higher rate. Also, note that the lifted product above has a much higher distance than the [[400,16,6]] hypergraph product constructed from the [16,4,6] classical LDPC code. 

# Bias-tailoring

## The CSS twisted toric code

Recall that the CSS toric code can be obtained 


```python
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



This code has periodic boundary conditions that link opposite sides of the lattice. The distance of this code can be increased by instead defing twisted boundary conditions. This can be achieved by taking the lifted product of two repetition code protographs:


```python
L=6
a1=pt.array([[(0,1)]])
a2=pt.array([[(0,2)]])

qcode=lifted_hgp(lift_parameter=L,a=a2,b=a1)
qcode.compute_code_distance() #this will take some time!
qcode.test()
```

    Warning: computing a code distance of codes with N>10 will take a long time.


    100%|██████████| 16383/16383 [00:00<00:00, 60090.39it/s]

    <Unnamed CSS code>, (2,4)-[[12,2,3]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -PCMs commute hx@hz.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unnamed CSS code> is a valid CSS code w/ params (2,4)-[[12,2,3]]


    





    True



Here we see that the boundary twist has increased the code distance from `d=2` to `d=3`.

### The CSS twisted toric code under infinite-bias

Now, consider the performance of the CSS twisted Toric code under X-bias. In this regime, the performance of the code depends exclusively on the `hz` component.


```python
h=qcode.hz
d=compute_code_distance(h)
print(f"Hz code distance = {d}")
```

    Hz code distance = 3


## The XZZX twisted toric code

The XZZX twisted toric code is obtained from the bias-tailored lifted product of two repetition codes. It is equivalent to the CSS twisted toric code up to a Hadamard rotation on the second block of `N/2` qubits. Using our code base, an XZZX twisted toric code can be constructed as follows


```python
from lifted_hgp import bias_tailored_lifted_product

L=6
a1=pt.array([[(0,1)]])
a2=pt.array([[(0,2)]])

qcode=bias_tailored_lifted_product(lift_parameter=L,a=a2,b=a1)
qcode.compute_code_distance() #this will take some time!
qcode.test()
```

    Warning: computing a code distance of codes with N>10 will take a long time.


    100%|██████████| 16383/16383 [00:00<00:00, 58158.03it/s]

    <Unamed stabiliser code>, [[12,2,3]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unamed stabiliser code> is a valid stabiliser code w/ params [[12,2,3]]


    





    True




```python
print(qcode.hx_proto)
```

    [[() (0,-2)]
     [(0,2) ()]]


### The XZZX twisted toric code under infinite bias
The XZZX toric code has improved distance in the infinite-bias regime. Eg. for the case of infinite X-bias, the `hz` code distance is


```python
h=qcode.hz
d=compute_code_distance(h)
print(f"Hz code distance = {d}")
```

    Hz code distance = 6


Not that this is an improvement on the infinite bias threshold of the CSS twisted toric code of the same size.

### General construction for XZZX twisted toric codes

In general, we define a twisted toric code on a `n1 x n2` lattice. The parity check matrices can be obtained via the bias-tailored lifted product. 


```python
from lifted_hgp import bias_tailored_lifted_product

n1=4
n2=3
L=n1*n2
a1=pt.array([[(0,1)]])
a2=pt.array([[(0,n2)]])

qcode=bias_tailored_lifted_product(lift_parameter=L,a=a2,b=a1)
# qcode.compute_code_distance() #this will take some time!
qcode.test()
```

    <Unamed stabiliser code>, [[24,2,nan]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unamed stabiliser code> is a valid stabiliser code w/ params [[24,2,nan]]





    True



## Bias-tailored LDPC codes

The bias-tailored lifted product can be used to create a quantum LDPC code from any pair of protographs. For example:


```python
a1=pt.array([
        [(0), (11), (7), (12)],
        [(1), (8), (1), (8)],
        [(11), (0), (4), (8)],
        [(6), (2), (4), (12)]])


qcode = bias_tailored_lifted_product(lift_parameter=13,a=a1,b=a1)
qcode.test()
```

    <Unamed stabiliser code>, [[416,18,nan]]
     -Block dimensions: Pass
     -PCMs commute hz@hx.T==0: Pass
     -lx \in ker{hz} AND lz \in ker{hx}: Pass
     -lx and lz anticommute: Pass
     -<Unamed stabiliser code> is a valid stabiliser code w/ params [[416,18,nan]]





    True



Now if we print the `hz` protograph, we see that it has been simplified to a set of 8 decoupled copies of original protograph `a1` (and its transpose) along the diagonal. Consequenlty, the quantum code in the infinite bias limit inherits the distance of the classical code `d=26`


```python
with np.printoptions(threshold=np.inf):
    print(qcode.hz_proto)
```

    [[(0) (11) (7) (12) () () () () () () () () () () () () () () () () () () () () () () () () () () () ()]
     [(1) (8) (1) (8) () () () () () () () () () () () () () () () () () () () () () () () () () () () ()]
     [(11) (0) (4) (8) () () () () () () () () () () () () () () () () () () () () () () () () () () () ()]
     [(6) (2) (4) (12) () () () () () () () () () () () () () () () () () () () () () () () () () () () ()]
     [() () () () (0) (11) (7) (12) () () () () () () () () () () () () () () () () () () () () () () () ()]
     [() () () () (1) (8) (1) (8) () () () () () () () () () () () () () () () () () () () () () () () ()]
     [() () () () (11) (0) (4) (8) () () () () () () () () () () () () () () () () () () () () () () () ()]
     [() () () () (6) (2) (4) (12) () () () () () () () () () () () () () () () () () () () () () () () ()]
     [() () () () () () () () (0) (11) (7) (12) () () () () () () () () () () () () () () () () () () () ()]
     [() () () () () () () () (1) (8) (1) (8) () () () () () () () () () () () () () () () () () () () ()]
     [() () () () () () () () (11) (0) (4) (8) () () () () () () () () () () () () () () () () () () () ()]
     [() () () () () () () () (6) (2) (4) (12) () () () () () () () () () () () () () () () () () () () ()]
     [() () () () () () () () () () () () (0) (11) (7) (12) () () () () () () () () () () () () () () () ()]
     [() () () () () () () () () () () () (1) (8) (1) (8) () () () () () () () () () () () () () () () ()]
     [() () () () () () () () () () () () (11) (0) (4) (8) () () () () () () () () () () () () () () () ()]
     [() () () () () () () () () () () () (6) (2) (4) (12) () () () () () () () () () () () () () () () ()]
     [() () () () () () () () () () () () () () () () (0) (-1) (-11) (-6) () () () () () () () () () () () ()]
     [() () () () () () () () () () () () () () () () (-11) (-8) (0) (-2) () () () () () () () () () () () ()]
     [() () () () () () () () () () () () () () () () (-7) (-1) (-4) (-4) () () () () () () () () () () () ()]
     [() () () () () () () () () () () () () () () () (-12) (-8) (-8) (-12) () () () () () () () () () () () ()]
     [() () () () () () () () () () () () () () () () () () () () (0) (-1) (-11) (-6) () () () () () () () ()]
     [() () () () () () () () () () () () () () () () () () () () (-11) (-8) (0) (-2) () () () () () () () ()]
     [() () () () () () () () () () () () () () () () () () () () (-7) (-1) (-4) (-4) () () () () () () () ()]
     [() () () () () () () () () () () () () () () () () () () () (-12) (-8) (-8) (-12) () () () () () () () ()]
     [() () () () () () () () () () () () () () () () () () () () () () () () (0) (-1) (-11) (-6) () () () ()]
     [() () () () () () () () () () () () () () () () () () () () () () () () (-11) (-8) (0) (-2) () () () ()]
     [() () () () () () () () () () () () () () () () () () () () () () () () (-7) (-1) (-4) (-4) () () () ()]
     [() () () () () () () () () () () () () () () () () () () () () () () () (-12) (-8) (-8) (-12) () () () ()]
     [() () () () () () () () () () () () () () () () () () () () () () () () () () () () (0) (-1) (-11) (-6)]
     [() () () () () () () () () () () () () () () () () () () () () () () () () () () () (-11) (-8) (0) (-2)]
     [() () () () () () () () () () () () () () () () () () () () () () () () () () () () (-7) (-1) (-4) (-4)]
     [() () () () () () () () () () () () () () () () () () () () () () () () () () () () (-12) (-8) (-8) (-12)]]


# BP+OSD decoding of bias-tailored LDPC codes


