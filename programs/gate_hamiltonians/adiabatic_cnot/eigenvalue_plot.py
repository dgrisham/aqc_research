#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt


def main():
    # define each piece of the total Hamiltonian
    Hi  = np.array( [[3,0,0,0], [0,2,0,0], [0,0,1,0], [0,0,0,0]], dtype=complex )
    Hf  = np.array( [[3,0,0,0], [0,2,0,0], [0,0,0,0], [0,0,0,1]], dtype=complex )
    Hif = np.array( [[0,0,0,0], [0,0,0,0], [0,0,0,1], [0,0,1,0]], dtype=complex )

    # total time to evolve through
    T = 1
    # time step
    step = 0.0001
    # arrays to hold results
    ## ground state energies
    evals0 = []
    ## first excited state energies
    evals1 = []
    ## energy gaps
    gap = []
    for t in np.arange(0, T, step):
        # get sorted list of eigenvalues
        evals = np.sort(np.linalg.eigvals(Htot(t/T, Hi, Hf, Hif)))
        # store ground state energy
        evals0 += [evals[0]]
        # store first excited state energy
        evals1 += [evals[1]]
        # calculate and store gap between ground state and first excited state
        gap += [evals[1] - evals[0]]

    # set y-axis plot range
    plt.ylim(0, 1)
    # time values for x-axis
    time_vals = np.arange(0, T, step)
    # plot stored values
    plt.plot(time_vals, evals0, color='blue')
    plt.plot(time_vals, evals1, color='purple')
    plt.plot(time_vals, gap, color='red')

    plt.show()

def Htot(s, Hi, Hf, Hif):
    return ((1-s)*Hi + s*Hf + s*(1-s)*Hif)

if __name__ == '__main__':
    main()
