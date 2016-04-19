#!/usr/bin/python3

# =============================================================================
# This script creates a local (3 qubit, 8x8) update operator, U, for the 16
# unitary ECAs and a unitary operator, V, to act on the center site if permitted
# by S.  I introduce a numbering scheme, S in [0, 1, ... 15]. The S relates to
# the original ECA code numer as (S, R):
#
# (0, 204), (1, 201), (2, 198),  (3, 195), (4, 156), (5, 153), (6, 150),
# (7, 147), (8, 108), (9, 105), (10, 102), (11, 99), (12, 60), (13, 157)
# (14, 57), (15, 51)
#
# Results are saved as .hdf5 and include single site reduced density matricies
# with key 'one_site', two site reduced density matrices : 'two_site'. The
# spin moments <Ai>, and <ABij> with A, B \in {X, Y, Z} are saved as symetric
# matrices. A != B has zero trace and A = B has <Ai> stored along the diagonal.
#
#
# By Logan Hillberry
# =============================================================================

import warnings
import copy
import time
import h5py

import numpy    as np
import simulation.fio      as io
import simulation.matrix   as mx
import simulation.states   as ss
import simulation.measures as ms

from os.path import isfile
from collections import OrderedDict
from math import fabs
from cmath import sqrt, pi, exp

# Constructing U
# ==============

# The center site gets updated with V if the neighbors are in a suitable
# configuration. Otherwise, the center site remains unchanged.


def make_V(V, s):
    V_conf = V.split('_')
    if len(V_conf) == 2:
        V_string, ph = V_conf
        ph = eval(ph)*pi/180.0
        Pmat = np.array( [[1.0,  0.0 ],[0.0 , exp(1.0j*ph)]], dtype=complex )
        ss.ops['P'] = Pmat 
    else:
        V_string = V_conf[0]
    Vmat= mx.listdot([ss.ops[k] for k in V_string])
    return s*Vmat + (1-s)*ss.ops['I']


def make_U(S, V, use_R=False, BC_conf='00'):
    # expand S into 4 digits of binary. 15 possible unitary ECAs
    # MSB comes first in the list: [s11, s10, s01, s00]
    Sb = np.array(mx.dec_to_bin(S, 4))

    if use_R is True:
        # first arg interpreted as the usual ECA code number

        # compute swap rule with XOR 204, expand and extract S
        R = np.array(mx.dec_to_bin(204^S, 8))
        S1= np.take(R, [0, 1, 4, 5])
        S2 = np.take(R, [2, 3, 6, 7])
        if not np.array_equal(S1, S2):
            # check if rule number is unitary
            warning_string = 'R = '+str(S)+' is not a unitary rule number'
            warnings.warn(warning_string)
            Sb = S1
        else:
            # set S
            Sb = S1

    # Order S into matrix [ [s00, s01], [s10, s11] ]
    S_mat = Sb[::-1].reshape((2,2)) # add .T for old R numbering (pre winter 2015)
    # prepare projectors for looping over in the order  {|0><0|, |1><1|}
    neighbor_projs = [ss.ops['0'], ss.ops['1']]

    # initialize U's
    U = np.zeros((8,8), dtype=complex)
    Ur = np.zeros((4,4), dtype=complex)
    Ul = np.zeros((4,4), dtype=complex)

    # loop through neighborhood configurations, applying V when logic of S permits
    for m, left_proj in enumerate(neighbor_projs):
        for n, right_proj in enumerate(neighbor_projs):
            Vmn = make_V(V, S_mat[m,n])
            U = U + mx.listkron([left_proj, Vmn, right_proj])

    # 2 qubit operator for right-most site, fix right boundary to |BC_conf[1]>
    for m, left_proj in enumerate(neighbor_projs):
        n=int(BC_conf[1])
        Vmn = make_V(V, S_mat[m,n])
        Ur = Ur + mx.listkron([left_proj, Vmn])

    # 2 qubit operator for left-most site, fix left boundary to |BC_conf[0]>
    for n, right_proj in enumerate(neighbor_projs):
        m=int(BC_conf[0])
        Vmn = make_V(V, S_mat[m,n])
        Ul = Ul + mx.listkron([Vmn, right_proj])

    return Ul, U, Ur


# Time evolution
# ==============

# update procedure for fixed BC's
# -------------------------------
def update_site(j, state, Ul, Uj, Ur, L, BC_typ = '1'):
    # site 0 has only one neighbor to the right
    if BC_typ is '1':
        if j == 0:
            js = [0,1]
            state = mx.op_on_state(Ul, js, state)

        # site L-1 has only one neighbor to the left
        elif j == L-1:
            js = [L-2, L-1]
            state = mx.op_on_state(Ur, js, state)

        # all other sites have left and right neighbors
        else:
            js = [(j-1), j, (j+1)]
            state = mx.op_on_state(Uj, js, state)

    elif BC_typ is '0':
        js = [(j-1)%L, j, (j+1)%L]
        state = mx.op_on_state(Uj, js, state)
    
    return state

# check normalization of state
# ----------------------------
def check_norm(state, t, tol):
    # calculate inner product
    ip = (state.conj().dot(state)).real

    # check if ip == 1 within supplied tolorance
    if fabs(ip - 1.0) < tol:
        return state

    else:
        warnings.warn('Non-normalized state at t = ' + str(t) + \
            ': <psi|psi> =' + str(ip) + ' has been normalized' )
        state = 1.0/sqrt(ip) * state
        return state


# construct generator for exact time evolved quantum state
# --------------------------------------------------------
def time_evolve(params, tol=1E-10, norm_check=False):
    # load simulation parameters
    L = params['L']
    T = params['T']
    mode = params['mode']
    BC = params['BC']
    V = params['V']

    if 'R' in params:
        use_R = True
        S = params['R']

    elif 'S' in params:
        use_R = False
        S = params['S']

    BC_typ = BC.split('_')[0]
    if BC_typ is '1':
        if len(BC.split('_')) == 1:
            BC_conf = '00'
        else:
            BC_conf = BC.split('_')[1]


    # make update operators for left/right boundaries and th bulk
    Ul, Uj, Ur = make_U(S, V, use_R=use_R, BC_conf=BC_conf)

    # If no state supplied, make from the IC param
    if 'init_state' in params:
        init_state = params['init_state']
    else:
        IC = params['IC']
        init_state = ss.make_state(L, IC)

    # yield the initial state
    state = np.array(init_state)
    yield state

    for t in range(T):

        # Sweep ECA
        if mode=='sweep':
            for j in range(L):
                state = update_site(j, state, Ul, Uj, Ur, L, BC_typ = BC_typ)

        # Alternating ECA (evens first)
        elif mode=='alt':
            for j in list(range(0, L, 2)) + list(range(1, L, 2)):
                state = update_site(j, state, Ul, Uj, Ur, L, BC_typ = BC_typ)

        # Block ECA
        elif mode=='block':
            for k in [0,1,2]:
                for j in range(k, L-1+k, 3):
                    if j!=L:
                        state = update_site(j, state, Ul, Uj, Ur, L, BC_typ = BC_typ)

        # don't check normalization by default
        if norm_check is True:
            state = check_norm(state, t, tol)
        # yield the updated state
        yield state


# make indices of sites in the smaller half of a bipartite cut
# ------------------------------------------------------------
def bi_partite_inds(L, cut):
    inds = [i for i in range(cut+1)]
    if cut > int(L/2):
        inds = list(np.setdiff1d(range(L), inds))
    return inds

# reduced density matrix of the smaller of all bi-partitions of the lattice
# -------------------------------------------------------------------------
def bi_partite_rdm(L, state):
    return np.array([mx.rdms(state, bi_partite_inds(L, cut))
            for cut in range(L-1)])

# compute the inverse participation ratio
# ---------------------------------------
def inv_participation_ratio(L, state):
    pr = np.sum(np.abs(state)**4)
    if pr == 0.0:
        return 0.0
    else:
        return 1.0 / pr

# import/create simulation results of one and two site reduced density matrices
# and all one and two point spin averages
# -----------------------------------------------------------------------------
def run_sim(params, force_rewrite = False,
        sim_tasks=['one_site', 'two_site', 'bi_partite', 'IPR']):

    if 'fname' in params:
        fname = params['fname']
    else:
        if 'IC' in params and params['IC'][0] == 'r':
            fname = io.make_file_name(params, iterate = False)
            #fname = io.make_file_name(params, iterate = True)
            # don't iterate file names with a unique IC name
        else:
            fname = io.make_file_name(params, iterate = False)
            # set the file name for each simulation
        params['fname'] = fname


    # see if sim has already ran
    f = h5py.File(fname, 'r')
    sim_ran = 'one_site' in f
    f.close()

    # check if file already exists and, if so, if it should be re-written
    if not sim_ran or force_rewrite:
        print('Running simulation...')

        # ensure one site rhos are calculated. That is how the sim_ran is determined
        if 'one_site' not in sim_tasks:
            sim_tasks.append('one_site')

        # collect params needed for initialization
        L = params['L']
        T = params['T']

        # initialize arrays for data collection
        data = {}
        if 'one_site' in sim_tasks:
            data['one_site'] = np.zeros((T+1, L, 2, 2), dtype = complex)

        if 'two_site' in sim_tasks:
            data['two_site'] = np.zeros((T+1, L, L, 4, 4), dtype = complex)

        if 'IPR' in sim_tasks:
            data['IPR'] = np.zeros(T+1)

        if 'bi_partite' in sim_tasks:
            # each cut is stored as a different data set of length T
            data['bi_partite'] = {}
            rdm_dims = [2**len(bi_partite_inds(L, cut))
                    for cut in range(L-1)]

            for cut, dim in enumerate(rdm_dims):
                data['bi_partite']['cut'+str(cut)] = np.zeros((T+1, dim, dim), dtype=complex)

        if 'sbond' in sim_tasks:
            data['sbond'] = np.zeros((T+1,L-1))
            data['snum'] = np.zeros((T+1,L-1))

        # loop through quantum states
        for t, state in enumerate(time_evolve(params)):
            # store the entire initial state
            if t == 0:
                data['init_state'] = state

            if 'IPR' in sim_tasks:
                data['IPR'][t] = inv_participation_ratio(L, state)

            # first loop through lattice, make single site matrices
            for j in range(L):
                if 'one_site' in sim_tasks:
                    rtj = mx.rdms(state, [j])
                    data['one_site'][t, j] = rtj

                # second loop through lattice
                if 'two_site' in sim_tasks:
                    for k in range(j+1, L):
                        rtjk = mx.rdms(state, [j, k])
                        data['two_site'][t, j, k] = rtjk
                        data['two_site'][t, k, j] = rtjk

            if 'bi_partite' in sim_tasks:
                for cut in range(L-1):
                    rtc = mx.rdms(state, bi_partite_inds(L, cut))
                    data['bi_partite']['cut'+str(cut)][t] = rtc

            if 'sbond' in sim_tasks:
                for b in range(L-1):
                    rtb = mx.rdms(state, bi_partite_inds(L, b))
                    sbond, snum = ms.vn_entropy(rtb, get_snum=True)
                    data['sbond'][t,b] = sbond
                    data['snum'][t,b] = snum 

        if 'IPR' in sim_tasks:
            freqs, amps, rn = ms.make_ft(data['IPR'])
            data['FIPR'] = amps
            data['RNIPR'] = rn
            data['freqs'] = freqs
    else:
        print('Importing states...')
        sim_tasks.append('FIPR')
        sim_tasks.append('freqs')
        data = io.read_hdf5(fname, sim_tasks)

    return data



# Default behavior of this file
# =============================
if __name__ == "__main__":
    import csv
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import simulation.plotting as pt
    import simulation.fitting as ft

    font = {'family':'serif', 'size':10}
    mpl.rc('font',**font)


    # Simulation time scaling with L
    # ------------------------------
    # set up loop to time 1 iteration of evolution for increasing L
    L_list = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
    t_list = []
    for L in L_list:

    # dicts of params specifying the simulation
        params =  {
                        'output_dir' : 'tmp/trash',

                        'L'    : L,
                        'T'    : 1,
                        'mode' : 'sweep',
                        'R'    : 102,
                        'V'    : ['H'],
                        'IC'   : 'l0',
                        'BC'   : '1'
                                           }


        # get run time of simulation, append to timing list
        tic = time.time()
        run_sim(params, force_rewrite = True, sim_tasks=['one_site', 'two_site', 'IPR'])
        toc = time.time()
        t_list.append(toc-tic)

    # save data to compare as improvements are made
    # NOTE: Change file names after each optimization!!
    data_fname = io.base_name('timing', 'data')+'nocut_iteration_timing2.csv'
    plots_fname = io.base_name('timing', 'plots')+'nocut_iteration_timing2.pdf'
    with open(data_fname, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(L_list, t_list))

    # check if data imports
    with open(data_fname, 'r') as f:
        reader = csv.reader(f)
        # csv reader imports numbers as strings, convert them to floats
        timing_list = np.asarray([[float(val) for val in row] \
                for row in list(reader)])

    # get data into two lists
    L_list = timing_list[:,0]
    t_list = timing_list[:,1]
    Ls = np.linspace(min(L_list), max(L_list), 150)


    # fit data to an exponential form
    Bs, chi2 = ft.f_fits(ft.fexp, [0.0, 2.0, 0.0], L_list, t_list)

    # plot the results
    plt.plot(L_list, t_list)
    plt.plot(Ls, ft.fexp(Bs, Ls),
            label = r"${%.6f} \times {%.2f}^L  {%+.3f}$" % (Bs[0], Bs[1], Bs[2])
            + "\n  $\chi^2 = {:.2e}$".format(chi2))
    plt.xlabel('L')
    plt.ylabel('time to simulate and save one iteration [s]')
    plt.grid('on')
    plt.legend(loc='upper left')
    plt.savefig(plots_fname)
    plt.close()
