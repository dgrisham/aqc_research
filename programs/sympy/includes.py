#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from sympy.physics.quantum.qubit import Qubit

q = Qubit('01101')
print("---includes---\n")
print(q)
print("len: {}\n".format(q.dimension))
print("--------------\n")
