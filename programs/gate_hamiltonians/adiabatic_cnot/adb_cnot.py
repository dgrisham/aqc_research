#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import qops.matrix as mx
import qops.states as st

sns.set_style("whitegrid")


def main():
    Hi  = np.array( [[3,0,0,0], [0,2,0,0], [0,0,1,0], [0,0,0,0]], dtype=complex )
    Hf  = np.array( [[3,0,0,0], [0,2,0,0], [0,0,0,0], [0,0,0,1]], dtype=complex )
    Hif = np.array( [[0,0,0,0], [0,0,0,0], [0,0,0,1], [0,0,1,0]], dtype=complex )

    state = np.array( [0,0,1,0], dtype=complex )

    T = 1
    step = 0.0001
    gap = []
    evals0 = []
    evals1 = []
    print("initial X measurement: {}".format(measure_all(state, 1, 'X')))
    print("initial Y measurement: {}".format(measure_all(state, 2, 'Y')))
    print("initial Z measurement: {}".format(measure_all(state, 2, 'Z')))
    for t in np.arange(0, T, step):
        state = mx.propagate(Htot(t/T, Hi, Hf, Hif), t, state)
        #print("X measurement: {}".format(measure_all(state, 2, 'X')))
        #print("Y measurement: {}".format(measure_all(state, 2, 'Y')))
        #print("Z measurement: {}".format(measure_all(state, 2, 'Z')))
        evals = np.sort(np.linalg.eigvals(Htot(t/T, Hi, Hf, Hif)))
        evals0 += [evals[0]]
        evals1 += [evals[1]]
        gap += [evals[1] - evals[0]]

    print("final X measurement: {}".format(measure_all(state, 2, 'X')))
    print("final Y measurement: {}".format(measure_all(state, 2, 'Y')))
    print("final Z measurement: {}".format(measure_all(state, 2, 'Z')))
    #gap = np.asarray(gap, dtype=float)

    #fig = plt.figure()
    plt.ylim(0, 1)
    plt.plot(np.arange(0, T, step), evals0, color='blue')
    plt.plot(np.arange(0, T, step), evals1, color='purple')
    plt.plot(np.arange(0, T, step), gap, color='red')
    plt.show()

###################################################################
## measure the spin along the z-axis for all qubits in the state ##
###################################################################
def measure_all(state, N, meas='Z', mode='rdms'):
    if mode == 'rdms':
        # reduced density matrix for each site calculated from the state
        rho_list = [mx.rdms(state, [j]) for j in range(N)]

        # measure z projection at each site. Take real part because measurements
        # of Hermitian ops always give a real result
        return [np.trace(rho.dot(st.pauli[meas])).real for rho in rho_list]
    else:
        results = []
        Z = st.pauli[meas]
        for i in range(N):
            results += [mx.measure(Z, i, state)]

        return results

def Htot(s, Hi, Hf, Hif):
    return ((1-s)*Hi + s*Hf + s*(1-s)*Hif)

def Htot_bad(s, Hi, Hf, Hif):
    return ((1-s)*Hi + s*Hf)

if __name__ == '__main__':
    main()
