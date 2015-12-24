#!/usr/bin/python

import numpy as np
from cmath import sqrt
from functools import reduce

# get a particular operator (choose 'I', 'X', 'Y', 'Z', or 'H')
def get_op(op):
    return get_ops()[op]

# returns a dictionary with all of our common operators
def get_ops():
    # initialize our general operators with the identity
    general_ops = identity()
    # add the pauli operators
    general_ops.update(paulis())
    # add the hadamard operator
    general_ops.update(hadamard())

    # return a dictionary of our operators, which are numpy arrays
    return { op : make_op(**args) for op, args in general_ops.items() }

# make a numpy array out of the input rows and overall multiplicative factor
def make_op(rows, factor=1):
    return factor * np.array([row for row in rows], dtype=complex)

# make a qubit with weights specified by cs, in the basis returned by
# basis()
def qubit(cs, factor=1):
    return factor * np.array([ci * bi for ci, bi in zip(cs, basis())])

# return the 2-dimensional computational basis states
def basis():
    basis = [[1,0], [0,1]]
    return(np.array(bi, dtype=complex) for bi in basis)

# return the commutator of A and B, i.e. [A,B]
def commutator(A, B):
    return dot(A, B) - dot(B, A)

# redefine the dot product so that it is called in the same way that
# sys_qubit_op is
def dot(op, q):
    # apply 'op' to q
    return op.dot(q)

# the cnot matrix -- will be used in get_opts (TODO: add this in)
def cnot():
    return {
        'cnot' : { 'rows' : ([1, 0, 0, 0],
                             [0, 1, 0, 0],
                             [0, 0, 0, 1],
                             [0, 0, 1, 0])
        }
    }

# the identity matrix -- used in get_ops()
def identity():
    return {
        'I' : { 'rows' :  ([1, 0],
                           [0, 1])
        }
    }

# the pauli matrices -- used in get_ops()
def paulis():
    return {
        'X' : { 'rows' :  ([0, 1],
                           [1, 0])
        },
        'Y' : { 'rows' :  ([0, -1j],
                           [1j, 0])
        },
        'Z' : { 'rows' :  ([1, 0],
                           [0, -1])
        }
    }

# the hadamard matrix -- used in get_ops()
def hadamard():
    return {
        'H' : { 'rows' :  ([1,1],
                           [1,-1]),
              'factor' :  1/sqrt(2)
        }
    }

# call main
if __name__ == '__main__':
    qs, ops = main()
