# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 22:25:10 2016

@author: user
"""

import MDAnalysis.analysis.distances as dist
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import ex_render

def run_rdf(tr, p_index):
    lx = 2
    ly = 2
    n_bins = 500
    rdfxy = np.zeros((n_bins, 2))
    frames_count = 1000
    tid = tr[0].particles.typeid == p_index
    n_par = len(tr[0].particles.position[tid, 0])
    p = np.zeros((n_par, 3), dtype=np.float32)
    for k in xrange(frames_count):
        frame = tr[k]
        if k % 10 == 0:
            print "Frame " + str(k) + " processing"
        p[:, :] = frame.particles.position[tid, 0:3]
        p[:, 2] = 0
        dists = dist.self_distance_array(p, box = np.array([lx, ly, 1, 90, 90, 90], dtype=np.float32))
        hist = np.histogram(dists, n_bins)
        rdfxy[:, 1] += hist[0]
    
    rdfxy[:, 1] /= 4*(2*np.pi*hist[1][0:-1]*(hist[1][1]-hist[1][0])*n_par/lx/ly)
    rdfxy[:, 1] *= 2.0/n_par/(n_par-1)
    rdfxy[:, 0] = hist[1][0:-1]-hist[1][0]
    return rdfxy

# Open trajectory file
file_name = "trajectory.gsd" #sys.argv[1]
tr = ex_render.get_trajectory(file_name)
## Select atoms
p_index = 1 #sys.argv[2]

rdf = run_rdf(tr, p_index)

#n_bins = 50
#r1 = np.zeros(n_bins-1)
#zmin = -5
#zmax = 5
#dz = 0.1
#zmas = np.arange(zmin, zmax, dz)
#u.atoms.translate([0, 0, zmin])
#for z in zmas:
#    print "z = " + str(z)
#    string = "resname DPPC and prop z > -"+str(dz)+" and prop z < "+str(dz)
#    sel = u.select_atoms(string)
#    r = rdf.InterRDF(sel, sel, nbins = n_bins, range=(0.0, 7.0))
#    r.run()
#    r1 += r.rdf[1:]
#    u.atoms.translate([0, 0, dz])
#
#plt.plot(r.bins[1:], r1/len(zmas)/n_bins)

#n_bins = 50
#string = "resname DPPC and name P"
#sel = u.select_atoms(string)
#r = rdf.InterRDF(sel, sel, nbins = n_bins, range=(0.0, 25.0))
#r.run()
#
#plt.plot(r.bins[1:], r.rdf[1:])