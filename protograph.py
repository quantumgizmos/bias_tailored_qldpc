import numpy as np
import copy as cp

def permutation_matrix(n: int,shift: int)->np.ndarray:
    '''
    Outputs a size-n permutation matrix.

    Parameters:
    -----------
        n: int
            matrix dimension
        shift: int
            the shift parameter
    
    Returns
    -------
        mat: nxn matrix shifted by `shift' columns to the left
    '''
    return np.roll(np.identity(n),1*shift,axis=1).astype(int)

class ring_of_circulants_f2():
    '''
    Class implementing the algebra of the ring of circulants over the field f2

    Parameters
    ----------
    non_zero_coefficients: int
        List of the non-zero terms in the polynomial expansion of the ring element
    '''
    
    def __init__(self,non_zero_coefficients):
        try:
            self.coefficients=list(non_zero_coefficients)
        except TypeError:
            self.coefficients=[non_zero_coefficients]
        self.coefficients=np.array(self.coefficients).astype(int)
        try:
            assert len(self.coefficients.shape)==1
        except AssertionError:
            raise TypeError("The input to ring_of_circulants_f2 must be a one-dimensional list")
        
    def __add__(self,x):
        '''
        Overload for the addition operator between two ring elements
        
        Parameters
        ----------
        self: ring_of_circulants_f2
        x: ring_of_circulants_f2

        Returns
        -------
        ring_of_circulants_f2
        '''
        return ring_of_circulants_f2(self.coefficients.tolist()+x.coefficients.tolist())
    
    def __repr__(self):
        return f"protograph.ring_of_circulants_f2({self.__str__()})"

    def __str__(self):
        '''
        What we see when we print()
        '''
        length=self.len()
        out="("
        for i,value in enumerate(self.coefficients):
            out+=str(value)
            if i != (length-1):
                out+=","
        out+=")"
        return out

    @property
    def T(self):
        '''
        Returns the transpose of an element from the ring of circulants

        Returns
        -------
        ring_of_circulants_f2
        '''
        transpose_coefficients=-1*self.coefficients
        return ring_of_circulants_f2(transpose_coefficients)

    def __mul__(self,other):
        '''
        Overloads the multiplication operator * between elements of the ring of circulants
        '''

        try:
            assert type(self)==type(other)
        except AssertionError:
            raise TypeError(f"Ring elements can only be multiplied by other ring elements. Not by {type(other)}")

        no_coeffs=self.len()*other.len()

        # print(no_coeffs)

        new_coefficients=np.zeros(no_coeffs).astype(int)
        for i,a in enumerate(self.coefficients):

            for j,b in enumerate(other.coefficients):
                new_coefficients[i*other.len() + j]=a+b
        
        return ring_of_circulants_f2(new_coefficients)

    def __len__(self):
        return len(self.coefficients)

    def len(self):
        return len(self.coefficients)
        

    def to_binary(self,lift_parameter):

        '''
        Converts ring element to its binary representation
        
        Parameters
        ----------
        lift_parameter:int
            The size of the permutation matrices used to map to binary
        
        Returns
        numpy.ndarray
            Binary matrix in numpy format
        '''

        mat=np.zeros((lift_parameter,lift_parameter)).astype(int)
        for coeff in self.coefficients:
            mat+=permutation_matrix(lift_parameter,coeff)
        return mat %2


class array():

    '''
    Class implementing a protograph (an array where the elements are in the ring of circulants)
    
    
    Parameters
    ----------
    proto_array: array_like, 2D
        The input should be of the form [[(0,1),(),(1)]] where each tuple is the input to the ring_of_circulants_f2 class
    '''

    def __init__(self,proto_array):

        # Reads in input arrays and converts tuples to ring_of_circulants_f2 objects
        temp_proto=np.array(proto_array).astype(object)
        if len(temp_proto.shape)==3:
            m,n,_=temp_proto.shape
        elif len(temp_proto.shape)==2:
            m,n=temp_proto.shape
        else:
            raise TypeError("The input protograph must be a three-dimensional array like object.")

        self.proto_array=np.empty((m,n)).astype(ring_of_circulants_f2)

        for i in range(m):
            for j in range(n):
                if isinstance(temp_proto[i,j],ring_of_circulants_f2):
                    self.proto_array[i,j]=temp_proto[i,j]
                else:
                    self.proto_array[i,j]=ring_of_circulants_f2(temp_proto[i,j])

    def __repr__(self):

        m,n=self.proto_array.shape
        out="[["

        for i in range(m):
            if i!=0:
                out+="["
            for j in range(n):
                out+=str(self.proto_array[i,j])
                if j!=n-1:
                    out+=","
            if i!=m-1:
                out+="],"
            else:
                out+="]]"

        return f"protograph.array({out})"

    def __str__(self):
        '''
        Generates what we see when we print
        '''
        
        m,n=self.proto_array.shape
        out="[["

        for i in range(m):
            if i!=0:
                out+=" ["
            for j in range(n):
                out+=str(self.proto_array[i,j])
                if j!=n-1:
                    out+=" "
            if i!=m-1:
                out+="]\n"
            else:
                out+="]]"

        return out

    def __getitem__(self,coords):
        '''
        Allows to access protograph elements using obj[i,j] 
        '''
        i,j=coords
        return self.proto_array[i,j]

    def __setitem__(self,coords,value):
        '''
        Allows us to assign protograph elements using obj[i,j]=(int,int)
        '''
        i,j=coords
        self.proto_array[i,j]=ring_of_circulants_f2(value)

    @property
    def T(self):
        '''
        Returns the transpose of the protograph
        '''
        m,n=self.proto_array.shape
        temp=np.copy(self.proto_array)
        for i in range(m):
            for j in range(n):
                temp[i,j]=temp[i,j].T
                
        return array(temp.T)

    @property
    def shape(self):
        '''
        Returns the shape of protograph
        '''
        return self.proto_array.shape

    def to_binary(self,lift_parameter):
        '''
        Converts the protograph to binary
        '''
        L=lift_parameter
        m,n=self.shape
        mat=np.zeros((m*L,n*L)).astype(int)
        for i in range(m):
            for j in range(n):
                mat[i*L:(i+1)*L,j*L:(j+1)*L]=self[i,j].to_binary(L)
        return mat


def hstack(proto_list):
    '''
    hstack funciton for protographs
    '''
    for i,proto in enumerate(proto_list):
        proto_list[i]=proto.proto_array
    temp=np.hstack(proto_list)
    return array(temp)

def vstack(proto_list):
    'vstack function for protographs'
    for i,proto in enumerate(proto_list):
        proto_list[i]=proto.proto_array
    temp=np.vstack(proto_list)
    return array(temp)


def zeros(size):
    '''
    Returns a protograph full of zero elements from the ring of circulants
    '''
    proto_array=np.zeros((size,size)).astype(object)
    for i in range(size):
        for j in range(size):
            proto_array[i,j]=np.array([])
    return array(proto_array)


def identity(size):
    '''
    Returns an identity protograph
    '''
    proto=zeros(size)
    for j in range(size):
        proto[j,j]={0}
    return proto


def kron(a,b):
    '''
    The tensor product between two protographs
    '''
    temp=np.kron(a.proto_array,b.proto_array)
    return array(temp)

def copy(a):
    '''
    Copies a protograph
    '''
    return cp.deepcopy(a)