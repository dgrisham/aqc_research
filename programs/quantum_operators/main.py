#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from lib.states import *
from lib.matrix import *

# example application of the functions in states.py and matrix.py
# author: David Grisham
#
# refer to README.md or README.pdf for more info

def main():
    # define input values #
    # L = number of qubits in our system
    L = 4
    # IC = initial condition for our system
    IC = 'qt180_p0'
    # give list of operators we want to apply to the qubits of interest
    OPs = [pauli['X'], pauli['Y']]
    # finally, give a list of the qubits we wish to apply OP to
    js = [1, 2]

    for i in range(1, L+1):
        sys = make_state(i, IC)
        print('state of sys w/ {} qubits: {}'.format(i, sys))

    init_state = make_state(L, IC)
    # process these inputs... #
    # make the initial state based on L, IC
    init_state = make_state(L, IC)
    # kronecker the OPs together to get a single operator
    OP = listkron(OPs)

    # apply the operator to our state #
    final_state = op_on_state(OP, js, init_state)

    # print out the L and IC...
    print("\nL = {L}, IC = {IC}\n".format(L=str(L), IC=str(IC)))
    # ...and print out the corresponding state that we just made
    print('initial state\n{}'.format(init_state))
    # and finally, 
    print('final state\n{}\n'.format(final_state))

if __name__ == '__main__':
    main()

