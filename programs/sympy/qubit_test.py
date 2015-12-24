#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sympy.physics.quantum.gate as gate
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.qapply import qapply

def main():
    return Qubit('00')1

def eqsup(n):
    state = Qubit('0'*n)

    for i in range(n):
        state = gate.H(i) * state

    return qapply(state)

if __name__ == '__main__':
    q = main()
