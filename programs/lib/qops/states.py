#!/usr/bin/python3


from cmath import sqrt, sin, cos, exp, pi
import numpy as np
import qops.matrix as mx


# Global constants
# ================
# dictionary of local operators and local basis (b for basis)
# -----------------------------------------------------------
# pauli spin matrices in the Z-basis
pauli = {
    'I' : np.array( [[1.0,  0.0 ],[0.0,   1.0]], dtype=complex ),
    'X' : np.array( [[0.0,  1.0 ],[1.0,   0.0]], dtype=complex ),
    'Y' : np.array( [[0.0, -1.0j],[1.0j,  0.0]], dtype=complex ),
    'Z' : np.array( [[1.0,  0.0 ],[0.0 , -1.0]], dtype=complex )
}
bvecs = {
    # computational basis state vectors
    '0'  : np.array( [1.0, 0.0], dtype=complex ),
    '1'  : np.array( [0.0, 1.0], dtype=complex ),
    # equal superposition
    'es' : np.array( [1./sqrt(2), 1./sqrt(2)], dtype=complex ),
}
# projectors onto the basis states
brhos = { 
    '0' : np.array( [[1.0,   0.0],[0.0,   0.0]], dtype=complex ),
    '1' : np.array( [[0.0,   0.0],[0.0,   1.0]], dtype=complex ),
}

ops = {}
ops.update(pauli)
ops.update(bvecs)


# Initial State Creation
# ======================
   
# Create Fock state
# -----------------
def fock(L, config, zero = '0', one = '1'):
    dec = int(config)
    state = [el.replace('0', zero).replace('1', one)
            for el in list('{0:0b}'.format(dec).rjust(L, '0'))]
    return mx.listkron([bvecs[key] for key in state])

# Create state with set of alive sites
# -------------------------------------
def alive_list(L, config):
    js = map(int, config.split('_'))
    js = [L - j - 1 for j in js]
    dec = sum((2**x for x in js))
    return fock(L, dec)

# Create state with one live sites
# --------------------------------
def one_alive(L, config):
    dec = 2**int(config)
    return fock(L, dec)

# Create states with equal superposition locally
# ----------------------------------------------
def es_list(L, config):
    js = map(int, config.split('_'))
    js = [L - j - 1 for j in js]
    dec = sum((2**x for x in js))
    return fock(L, dec, one='es')

# Create state with all sites living
# ----------------------------------
def all_alive(L, config):
    dec = sum ([2**n for n in range(0,L)])
    return fock(L, dec)

# Create GHZ state
# ----------------
def ghz(L, congif):
    s1=['1']*(L)
    s2=['0']*(L)
    return (1.0/sqrt(2.0)) \
            * ((mx.listkron([bvecs[key] for key in s1]) \
                + mx.listkron([bvecs[key] for key in s2])))

# Create W state
# --------------
def W(L, config):
    return (1.0/sqrt(L)) \
            * sum ([one_alive(L, k) for k in range(L)])

# Create a state with GHZ-type entanglement.
# Reduces to 1/sqrt(2) (|00> + |11>) in L = 2 limit
# -------------------------------------------------
def entangled_list(L, config):
    js = map(int, config.split('_'))
    js = [L - j - 1 for j in js]
    dec = sum((2**x for x in js))
    return 1./sqrt(2) * (fock(L, 0) + fock(L, dec))

# Create a state with random single-qubit states
# ----------------------------------------------
def rand_state(L, config):
    # default background state (DOWN state)
    bg = np.array([0, 1], dtype=complex)
    # default state to be distribute among the background
    state = np.array([1, 0], dtype=complex)
    # get the probability specified in config
    state_prob = eval('.' + config)
    # array to hold the probabilities of state and bg
    prob = [state_prob, 1 - state_prob]
    # get random distribution of states based on prob
    states = [state, bg]
    distrib = np.random.choice([0, 1], size=L, p=prob)

    # kronecker the states in the distribution and return
    return mx.listkron([states[i] for i in distrib])

# Create a state with any of the above states 
# embeded in the center of the lattice
# -------------------------------------------
def center(L, config):
    len_cent = int(config[0])
    len_back = L - len_cent
    len_L = int(len_back/2)
    len_R = len_back - len_L
    cent_IC = [(config[1:], 1.0)]
    left = fock(len_L, 0)
    cent = make_state(len_cent, cent_IC)
    right = fock(len_R, 0)
    if len_back == 0:
        return cent
    elif len_back == 1:
        return mx.listkron([cent, right])
    else:
        return mx.listkron([left, cent, right])

def qubit(th, ph):
    return cos(th/2.0) * bvecs['0'] + exp(1j*ph) * sin(th/2) * bvecs['1']

def qubits(L, config):
    Tt, Pp = config.split('_')
    ang_dict = {
        'T' : np.linspace(0.0,  pi*float(Tt[1:]), L),
        't' : [float(Tt[1:])*pi/180.0]*L,
        'P' : np.linspace(0.0, 2*pi*float(Pp[1:]), L),
        'p' : [float(Pp[1:])*pi/180.0]*L,
    }
    th_list = ang_dict[Tt[0]]
    ph_list = ang_dict[Pp[0]]
    qubit_list = [0.0]*L
    for j, (th, ph) in enumerate(zip(th_list, ph_list)):
        qubit_list[j] = qubit(th, ph)
    return mx.listkron(qubit_list)

# Make the specified state
# ------------------------
smap = { 'd' : fock,
         's' : es_list,
         'l' : alive_list,
         'a' : all_alive,
         'c' : center,
         'q' : qubits,
         'G' : ghz,
         'W' : W,
         'E' : entangled_list,
         'r' : rand_state
}

def make_state(L, IC):
    state = np.array([0.0]*(2**L), dtype = complex)

    if type(IC) == str:
        name = IC[0]
        config = IC[1:]
        state = smap[name](L, config)

    if type(IC) == list:
        for s in IC: 
            name = s[0][0]
            config = s[0][1:]
            coeff = s[1]
            state = state + coeff * smap[name](L, config)

    state = mx.edit_small_vals(state)
    return state

def magnitude(state):
    return sum(np.absolute(state[i]) for i in range(len(state)))


def normalize(state):
    return (state/magnitude(state))

