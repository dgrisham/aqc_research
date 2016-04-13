#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

# author: David Grisham
# email: dgrisham@mines.edu
#
# requires: python3, custom qops libraries (included), numpy, scipy
#
# program to find a solution to the number partitioning problem using quantum
# adiabatic optimization. the algorithm here draws from:
#   - 'Ising Formulations of many NP problems' by Andrew Lucas,
#   - 'Solving NP-Hard Problems on an Adiabtic Quantum Computer' by Dan Padilha


import numpy as np
import qops.matrix as mx
import qops.states as st
from random import random
from time import sleep


###################################################################################
## main function; will run when this program is called from a shell/command line ## 
###################################################################################
def main():
    # set of values we would like to partition
    vals = [2, 4, 9, 15]

    print("\nthe values to be partitioned: {}\n".format(vals))

    # total time to evolve over
    T = 10
    # time step
    dt = 0.1
    # perform the adiabatic evolution
    final_state, info = evolve(vals, T, dt)
    # measure the spins of each qubit in the resulting state
    results = measure_all(final_state, len(vals) - 1)

    print("\nmeasured spins after adiabatic evolution: {}".format(results))

    return J, h, vals, final_state, results, info

def evolve(vals, T, dt):
    # strip off the first value as the constant field term
    field = vals[0]
    vals = vals[1:]
    # get the number of spins
    N = len(vals)

    print("\nthe value {} will be used as the constant field term, while each of the\n"\
            "other values ({}) will correspond to a qubit in our system\n".
        format(field, vals))

    # create the initial hamiltonian
    Hinit = build_Hinit(N)
    # create the final hamiltonian to evolve to
    Hfinal = build_Hfinal(field, vals, N)

    # make the initial state as the ground state of Hinit
    IC = 'qt90_p0'
    state = st.make_state(N, IC)

    # create a dictionary to hold some information for us
    info = {
        'states' : [state],
        'measurements' : [measure_all(state, N)],
        'Hinit' : Hinit,
        'Hfinal': Hfinal
    }

    print("final hamiltonian:\n{}\n".format(Hfinal))

    # apply the hamiltonian for each time step
    for t in np.arange(T, step=dt):
        # calculate the hamiltonian based on the current time
        H = (1 - t/T) * Hinit + (t/T) * Hfinal
        # apply the propagator to the state
        state = mx.propagate(H, dt, state)

        ## these lines do a bit of book-keeping for us
        # save calculated state in info
        info['states'] += [state]
        # measure the current state and save in info as well
        info['measurements'] += [measure_all(state, N)]

    return state, info

###################################
## build the initial hamiltonian ##
###################################
def build_Hinit(N):
    # define the pauli X matrix as 'X' for convenience
    X = st.pauli['X']
    # form Hinit (see Lucas paper for more info)
    Hinit = np.zeros((2**N, 2**N), dtype=complex)
    for i in range(N):
        Hinit += mx.make_big_mat([X], [i], N)

    return -1*Hinit

#################################
## build the final hamiltonian ##
#################################
def build_Hfinal(field, vals, N):
    return -1*(build_H_field(field, vals, N) + build_H_interaction(vals, N))

###################################################
## build the field term of the final hamiltonian ##
###################################################
def build_H_field(field, vals, N):
    # define the pauli Z matrix as 'Z' for convenience
    Z = st.pauli['Z']
    # create an empty Hfinal_h that we can fill in
    Hfinal_h = np.zeros((2**N, 2**N), dtype=complex)

    print("terms in the 'h' matrix:")
    for i in range(N):
        # get the h value for the current i
        hi = h(field, vals, i)
        print("h[{}]: {}".format(i, hi))
        # add in corresponding term to Hfinal_h
        Hfinal_h += mx.make_big_mat(hi * [Z], [i], N)

    return Hfinal_h

#########################################################
## build the interaction term of the final hamiltonian ##
#########################################################
def build_H_interaction(vals, N):
    # define the pauli Z matrix as 'Z' for convenience
    Z = st.pauli['Z']
    # kronecker two pauli Z's to operate on qubit pairs
    OP = mx.listkron([Z, Z])

    # create an empty Hfinal_J that we can fill in
    Hfinal_J = np.zeros((2**N, 2**N), dtype=complex)

    print("\nnon-zero terms in the 'J' matrix:")
    for j in range(N):
        for i in range(j):
            # get the J value for the current i,j pair
            Jij = J(vals, i, j)
            print("J[{}][{}]: {}".format(i,j, Jij))
            # add in corresponding term to Hfinal_J
            Hfinal_J += mx.make_big_mat(Jij * [Z, Z], [i, j], N)

    return Hfinal_J


############################################################
## get the [i,j]th entry of J, an upper-triangular matrix ##
############################################################
def J(vals, i, j):
    N = len(vals)
    if not ((0 <= i < (N - 1)) and ((i + 1) <= j < N)):
        # if i,j are not in the correct range, return 0
        return 0
    else:
        # else, return the corresponding entry in J
        return 2 * vals[i] * vals[j]

##############################################
## get the i'th value of h, a column vector ##
##############################################
def h(field, vals, i):
    return 2 * field * vals[i]


###################################################################
## measure the spin along the z-axis for all qubits in the state ##
###################################################################
def measure_all(state, N, mode='rdms'):
    if mode == 'rdms':
        # reduced density matrix for each site calculated from the state
        rho_list = [mx.rdms(state, [j]) for j in range(N)]

        # measure z projection at each site. Take real part because measurements
        # of Hermitian ops always give a real result
        return [np.trace(rho.dot(st.pauli['Z'])).real for rho in rho_list]
    else:
        results = []
        Z = st.pauli['Z']
        for i in range(N):
            results += [mx.measure(Z, i, state)]

        return results


####################################################################
## calculate the expected energy of the system (currently unused) ##
####################################################################
def expected_energy(vals):
    return sum(val**2 for val in vals)*-1


if __name__ == '__main__':
    J, h, vals, final_state, results, info = main()

