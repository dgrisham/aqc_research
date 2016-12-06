#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

#   author: David Grisham
#   email: dgrisham@mines.edu
#
#   requires: python3, custom qca libraries (included in repo), numpy, scipy,
#   matplotlib
#
#   program to find a solution to the number partitioning problem using quantum
#   adiabatic optimization. the algorithm here draws from:
#     - 'Ising Formulations of many NP problems' by Andrew Lucas,
#     - 'Solving NP-Hard Problems on an Adiabtic Quantum Computer' by Dan Padilha
#
#   to run:
#       `./partitioning.py "i1 i2 i3 ... in`
#   where `i1, i2, ..., in` are the numbers to partition
#   for example, to partition [1,4,9,15], we would run:
#       `./partitioning.py "1 4 9 15"`
#


import sys
import numpy as np
import matplotlib.pyplot as plt
import qca.simulation.matrix as mx
import qca.simulation.states as st
from random import random
from time import sleep

DEBUG = 1
DEBUG_MORE = 0
DEBUG_EVEN_MORE = 0


###################################################################################
## main function; will run when this program is called from a shell/command line ## 
###################################################################################
def main(argv, strip_field=True):

    # check for user input
    if len(argv) != 2 or argv[1] in ('-h', '--help'):
        print("usage: ./partitioning.py 'i1 i2 ... in'")
        exit()

    # set of values we would like to partition
    vals = [int(i) for i in argv[1].split()]
    print("input vals: {}".format(vals))

    print("\n***************************************")

    print("\nthe values to be partitioned: {}\n".format(vals))

    if strip_field:
        # strip off the first value as the constant field term
        field = vals[0]
        vals  = vals[1:]

        print("\nthe value {} will be used as the constant field term,"\
            "while each of the\n other values ({}) will correspond to a"\
            "qubit in our system\n".format(field, vals))
    else:
        field = 0
        print("\nno field term stripped for this run (all h_i = 0)")

    # total time to evolve over
    T = 1
    # time step
    dt = 0.001
    # perform the adiabatic evolution
    final_state, info = evolve(vals, field, T, dt, strip_field)
    # measure the spins of each qubit in the resulting state
    final_measurement = measure_all(final_state, len(vals))
    #final_measurement = measure_all((1/np.sqrt(2))*np.array([1,1]), 1)

    print("\nmeasured spins after adiabatic evolution: {}".
            format(final_measurement))
    print("\n***************************************")

    return J, h, vals, final_state, final_measurement, info

###############################################################
## built the Hamiltonian and adiabatically evolve the system ## 
###############################################################
def evolve(vals, field, T, dt, strip_field=True):
    if DEBUG:
        if strip_field: print_h(field, vals)
        print_J(vals)

    # get the number of spins
    N = len(vals)

    # create the initial hamiltonian
    Hinit = build_Hinit(N)
    #Hinit = Hinit / np.linalg.norm(Hinit)
    print("Hinit as returned:\n{}".format(Hinit))
    # create the final hamiltonian to evolve to
    Hfinal = build_Hfinal(field, vals, N, strip_field)
    #Hfinal = Hfinal / np.linalg.norm(Hfinal)

    # make the initial state as the ground state of Hinit
    IC = 'st90-p180'
    state = st.make_state(N, IC)
    #state = get_ground_state(Hinit, N)[1]
    print("initial state:\n{}".format(state))

    # get Hinit's energy, ground state vector, and measurement
    init_gs_energy, init_gs_vector, init_gs_measurement =\
        get_ground_state(Hinit, N)

    # create a dictionary to hold some information for us
    info = {
        'states' : [state],
        'measurements' : [measure_all(state, N)],
        'energies' : [expectation_val(state, N, Hinit)],
        'gs_energies' : [init_gs_energy],
        'gs_vectors' : [init_gs_vector],
        'gs_measurements' : [init_gs_measurement],
        'Hinit' : Hinit,
        'Hfinal': Hfinal,
        'Htot'  : [Hinit]
    }

    if DEBUG_MORE:
        print("\nfinal hamiltonian:\n{}\n".format(Hfinal))

    # apply the hamiltonian for each time step
    for t in np.arange(T, step=dt):
        # calculate the hamiltonian based on the current time
        if ((T-t) < dt):
            H = Hfinal
        else:
            H = (1 - t/T) * Hinit + (t/T) * Hfinal

        info['Htot'] += [H]
        # apply the propagator to the state
        state = mx.propagate(H, dt, state)

        ## these lines do a bit of book-keeping for us
        # save calculated state in info
        info['states'] += [state]
        # measure the current state and save in info as well
        info['measurements'] += [measure_all(state, N)]
        # store the expectation value of <state | H | state>
        info['energies'] += [expectation_val(state, N, H)]
        # save the groundstate energy and vector for the current H(s)
        gs_energy, gs_vector, gs_measurement = get_ground_state(H, N)
        info['gs_energies'] += [gs_energy]
        info['gs_vectors'] += [gs_vector]
        info['gs_measurements'] += [gs_measurement]

    plot_eigvals(info, T, dt, strip_field, save=False)

    return state, info

def expectation_val(state, N, operator):
    return (np.conj(state)).dot(operator).dot(state)

def get_ground_state(H, N):
    # get the eigenvalues
    eig_vals, eig_vecs = np.linalg.eig(H)
    # get the indices of the eigenenergies when sorted
    eig_order = np.argsort(eig_vals)
    # get ground state energy value
    gs_energy = eig_vals[eig_order[0]]
    # get ground state vector
    gs_vector = eig_vecs[eig_order[0]]
    # measure the qubits of this state
    gs_measurement = measure_all(gs_vector, N)

    return gs_energy, gs_vector, gs_measurement


def plot_eigvals(info, T, dt, strip_field=True, save=False):
    # arrays to hold results
    ## ground state energies
    evals0 = []
    ## first excited state energies
    evals1 = []
    ## energy gaps
    gap = []

    state_energies = []
    for t in np.arange(0, T, step=dt):
        # get sorted list of eigenvalues
        evals = np.sort(np.linalg.eigvals(info['Htot'][int(t/dt)]))
        # store ground state energy
        evals0 += [evals[0]]
        # store first excited state energy
        evals1 += [evals[1]]
        # calculate and store gap between ground state and first excited state
        gap += [np.abs(evals[1] - evals[0])]

        state_energies += [info['energies'][int(t/dt)]]

    # set y-axis plot range
    #plt.ylim(0, 1)
    # time values for x-axis
    time_vals = np.arange(0, T, dt)

    # plot stored values
    if strip_field:
        plt.title("Energy vs. Time (h_i field term set)")
    else: 
        plt.title("Energy vs. Time (only Jij terms, no h_i)")

    plt.plot(time_vals, evals0, color='blue')
    #plt.plot(time_vals, evals1, color='purple')
    plt.plot(time_vals, state_energies, color='green')
    #plt.plot(time_vals, gap, color='red')

    if save:
        plt.savefig("npp_adiabatic_energy-v-time.png")
    else:
        plt.show()

###################################
## build the initial hamiltonian ##
###################################
def build_Hinit(N):
    # define the ops X matrix as 'X' for convenience
    X = st.ops['X']
    # create Hinit (see Lucas paper for more info)
    Hinit = np.zeros((2**N, 2**N), dtype=complex)
    for i in range(N):
        Hinit += mx.make_big_mat([X], [i], N)

    h0 = 4500 # works for [2,4,9,15]
    #h0 = 5 # works for [1,2,3]
    return h0 * Hinit

#################################
## build the final hamiltonian ##
#################################
def build_Hfinal(field, vals, N, strip_field=True):
    H_interaction = build_H_interaction(vals, N)

    if strip_field:
        H_field = build_H_field(field, vals, N)

        if DEBUG_MORE:
            print("Hfield: {}".format(H_field))
            print("Hinteraction: {}".format(H_interaction))

        return -1 * (H_interaction + H_field)
    else:
        # no field term, don't add H_field to final Hamiltonian

        if DEBUG_MORE:
            print("Hinteraction: {}".format(H_interaction))

        return -1 * H_interaction

###################################################
## build the field term of the final hamiltonian ##
###################################################
def build_H_field(field, vals, N):
    # define the ops Z matrix as 'Z' for convenience
    X = st.ops['X']
    # create an empty Hfinal_h that we can fill in
    Hfinal_h = np.zeros((2**N, 2**N), dtype=complex)

    print("terms in the 'h' matrix:")
    for i in range(N):
        # get the h value for the current i
        hi = h(field, vals, i)
        Xi = hi*X
        #print("h[{}]: {}".format(i, Xi))
        # add in corresponding term to Hfinal_h
        Hfinal_h += mx.make_big_mat([Xi], [i], N)

    return Hfinal_h

#########################################################
## build the interaction term of the final hamiltonian ##
#########################################################
def build_H_interaction(vals, N):
    # define the ops Z matrix as 'Z' for convenience
    Z = st.ops['Z']

    # create an empty Hfinal_J that we can fill in
    Hfinal_J = np.zeros((2**N, 2**N), dtype=complex)

    for j in range(N):
        for i in range(j):
            # get the J value for the current i,j pair
            Jij = J(vals, i, j)
            Zi = Zj = Jij*Z
            # add in corresponding term to Hfinal_J
            if DEBUG_EVEN_MORE:
                print("\nadding to Hfinal_J w/ Jij ="
                        "{}:\n{}\n".format(Jij, Zi))
            Hfinal_J += mx.make_big_mat([Zi, Zj], [i, j], N)

    return Hfinal_J


############################################################
## get the [i,j]th entry of J, an upper-triangular matrix ##
############################################################
def J(vals, i, j):
    N = len(vals)
    if not ((0 <= i < (N - 1)) and ((i+1) <= j < N)):
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

########################
## print the J matrix ##
########################
def print_J(vals):
    N = len(vals)
    Jmat = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            Jmat[i][j] = J(vals, i, j)

    print("\nJ-matrix:\n{}".format(Jmat))

########################
## print the h vector ##
########################
def print_h(field, vals):
    N = len(vals)
    hvec = np.zeros(N)

    for i in range(N):
        hvec[i] = h(field, vals, i)

    print("\nh-vector\n{}".format(hvec))

####################################################################
## measure the spin along the z-axis for all qubits in the system ##
####################################################################
def measure_all(state, N, op=st.ops['Z']):
    # reduced density matrix for each site calculated from the state
    rho_list = [mx.rdms(state, [j]) for j in range(N)]
    #print("\nrho trace: {}".format(np.trace(rho)))

    # measure z projection at each site. Take real part because measurements
    # of Hermitian ops always give a real result
    return [np.trace(rho.dot(op)).real for rho in rho_list]

####################################################################
## calculate the expected energy of the system (currently unused) ##
####################################################################
def expected_energy(vals):
    return sum(val**2 for val in vals)*-1


if __name__ == '__main__':
    # run with user input
    J, h, vals, final_state, final_measurement, info = main(sys.argv, strip_field=False)

