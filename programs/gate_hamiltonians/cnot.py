#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from itertools import combinations


def main():
    n = 4 # number of qubits
    # define the characteristic cnot states
    cnot_states = ['0000', '0101', '1011', '1110']
    # break states down into equations of the form sum_i(hi) + sum_ij(Jij)
    eqs = get_eqs(cnot_states)

    print("states:\n{}".format(eqs))

def get_eqs(states):
    # array to hold our equations
    eqs = np.zeros((4, 10))

    # get equation for each state
    for i, state in enumerate(states):
        eqs[i] = get_eq(state)

    return eqs

def get_eq(state):
    # array to hold our equation
    eq = np.zeros(10)
    # dictionary to give us a mapping from state to eigenvalue of sigz
    eigenvals = {
            '0' : 1,
            '1' : -1
    }
    # get the signs for the hi's in the equation
    for i, q in enumerate(state):
        eq[i] = eigenvals[q]
    # get the signs for the Jij's in the equation
    for i, pair in enumerate(combinations(state, 2)):
        eq[i+4] = get_eigen(pair, eigenvals)

    return eq

def get_eigen(pair, eigenvals):
    result = 1
    # multiply the eigenvals for each qubit together to get the eigenvalue of
    # the state
    for q in pair:
        result *= eigenvals[q]

    return result

if __name__ == '__main__':
    main()

