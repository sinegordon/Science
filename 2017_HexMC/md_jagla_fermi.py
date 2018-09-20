#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 21:17:45 2017

@author: user
"""

import numpy as np
import hoomd, hoomd.md


def lj(r, rmin, rmax, epsilon, sigma):
    V = 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6);
    F = 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
    return (V, F)

def jf(r, rmin, rmax, n, epsilon, a, A0, A1, A2, B0, B1, B2):
    V = epsilon * ( (a / r)**n + A0/(1 + np.exp(A1/A0*(r/a-A2))) - B0/(1 + np.exp(B1/B0*(r/a-B2))))
    F = epsilon * ( -12*(a / r)**12/r - A1*np.exp(A1/A0*(r/a-A2))/a/(1+np.exp(A1/A0*(r/a-A2)))**2 +
                   B1*np.exp(B1/B0*(r/a-B2))/a/(1+np.exp(B1/B0*(r/a-B2)))**2)
    return (V, F)



hoomd.context.initialize()
unitcell=hoomd.lattice.sc(a=2.0, type_name='A')
hoomd.init.create_lattice(unitcell=unitcell, n=10)
nl = hoomd.md.nlist.cell()

table = hoomd.md.pair.table(width=1000, nlist=nl)
table.pair_coeff.set('A', 'A', func=jf, rmin=0.8, rmax=3.0, 
                     coeff=dict(n=20, epsilon=1.0, a=1.0, A0=4.56, A1=28.88, A2=1.36, B0=1.0, B1=3.57, B2=2.36))

all = hoomd.group.all();
d = hoomd.dump.gsd("md.gsd", period=10, group=all, overwrite=True)
hoomd.md.integrate.mode_standard(dt=0.005)
hoomd.md.integrate.langevin(group=all, kT=1.2, seed=4)
hoomd.run(10e3)