# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 22:25:10 2016

@author: user
"""

import MDAnalysis.analysis.distances as dist
import matplotlib.pyplot as plt
import numpy as np
import time
import ex_render
import sys


# Open trajectory file
file_name = sys.argv[1]
tr = ex_render.get_trajectory(file_name)
## Select atoms
p_index = sys.argv[2]


n_bins = 100
rdf = np.zeros(n_bins)
k = 0
frames_count = 100
for frame in tr:
    print "Processing frame # ", k
    #v = frame.velocities[inds]
    tid = frame.particles.typeid == p_index
    p = frame.particles.position[tid,0:2]
    dists = dist.self_distance_array(p, backend='OpenMP')
    hist = np.histogram(dists, n_bins)
    rdf += hist[0]
    k += 1
n_par = len(p)
rdf /= frames_count
rdf /= ( n_par / 2)
rdf /= (4*np.pi*hist[1][0:-1]**2*(hist[1][1]-hist[1][0])*n_par/u.dimensions[0]/u.dimensions[1]/u.dimensions[2])
plt.plot(hist[1][0:100], rdf[0:100])
#plt.plot(r.bins[1:], r.rdf[1:])

# Select tail
#string_tail = "(name C2* or name C3*) and not name C2 and not name C3 and not name C21 and not name C31"
#string_N = "name N"
#string_P = "name P"
#strings = {string_N : "N", string_P: "P", string_tail: "tail"}


#for frame in u.trajectory:
#    print sel.indices

# Select indices of integestng atoms
#string0 = "name C9 "
#for n in range(13, 41):
#    string0 = string0 + "or name C" + str(n) + " "

#n_bins = 50
#rdfxy = np.zeros(n_bins)
#zmin = -5
#zmax = 5
#dz = 0.1
#zmas = np.arange(zmin, zmax, dz)
#for z in zmas:
#    print "z = " + str(z)
#    for frame in u.trajectory:
#        string = "prop z > " + str(z) + " and prop z < " + str(z+dz)
#        sel = u.select_atoms(string)
#        inds = sel.indices
#        #v = frame.velocities[inds]
#        p = frame.positions[inds]
#        #inds = np.where(np.abs(p[:,2])<z)
#        p[:, 2] = 0
#        dists = dist.self_distance_array(p, box = u.dimensions)
#        hist = np.histogram(dists, n_bins)
#        rdfxy += hist[0]/(hist[1][0:-1])
#
##rdfxy[0] = 0
#plt.plot(hist[1][0:49], 2*rdfxy[0:49]/len(zmas)/len(p)/len(p))
#def run_rdf(temp, strings, string, mode):
#    if mode == 'all':
#        # Open trajectory file
#        bt = time.time()
#        TPR = './'+str(temp)+'/charmm-gui/gromacs/step7_1.tpr'
#        TRR = './'+str(temp)+'/charmm-gui/gromacs/step7_1nj.trr'
#        u = ma.Universe(TPR, TRR)
#        ma.core.flags['use_periodic_selections'] = True
#        ma.core.flags['use_pbc'] = False
#        et = time.time()
#        print "Open timing is ", et-bt
#        ## Select atoms
#        string = "resname DPPC and ( " + string + " )"
#        sel = u.select_atoms(string)
#        inds = sel.indices
#        n_bins = 100
#        rdfxy = np.zeros((n_bins, 2))
#        frames_count = 10000
#        n_par = len(inds)
#        p = np.zeros((n_par, 3), dtype=np.float32)
#        #rsd = [u.select_atoms("resid "+str(i) + " and (" + string + ")") for i in xrange(0, n_lipids)]
#        for k in xrange(frames_count):
#            frame = u.trajectory[k]
#            if k % 10 == 0:
#                print "Frame # " + str(k) + " processing"
#            #for i in xrange(0, n_lipids):
#            #    p[i, :] = rsd[i].center_of_mass()[:]
#            p[:, :] = frame.positions[inds, :]
#            p[:, 2] = 0
#            dists = dist.self_distance_array(p, box = u.dimensions)
#            hist = np.histogram(dists, n_bins)
#            rdfxy[:, 1] += hist[0]
#        
#        rdfxy[:, 1] /= 4*(2*np.pi*hist[1][0:-1]*(hist[1][1]-hist[1][0])*n_par/u.dimensions[0]/u.dimensions[1])
#        rdfxy[:, 1] *= 2.0/n_par/(n_par-1)
#        rdfxy[:, 0] = hist[1][0:-1]-hist[1][0]
#        return rdfxy
#        
#for temp in [298, 310, 318, 333]:
#    for s in strings:
#        print s
#        rdf = run_rdf(temp, strings, s, 'all')
#        np.savetxt("./RDF/rdfxy_"+str(temp)+"_"+strings[s]+".xvg", rdf)
#        #plt.plot(rdf[:, 0], rdf[:, 1])
#
##n_bins = 50
##r1 = np.zeros(n_bins-1)
##zmin = -5
##zmax = 5
##dz = 0.1
##zmas = np.arange(zmin, zmax, dz)
##u.atoms.translate([0, 0, zmin])
##for z in zmas:
##    print "z = " + str(z)
##    string = "resname DPPC and prop z > -"+str(dz)+" and prop z < "+str(dz)
##    sel = u.select_atoms(string)
##    r = rdf.InterRDF(sel, sel, nbins = n_bins, range=(0.0, 7.0))
##    r.run()
##    r1 += r.rdf[1:]
##    u.atoms.translate([0, 0, dz])
##
##plt.plot(r.bins[1:], r1/len(zmas)/n_bins)
#
##n_bins = 50
##string = "resname DPPC and name P"
##sel = u.select_atoms(string)
##r = rdf.InterRDF(sel, sel, nbins = n_bins, range=(0.0, 25.0))
##r.run()
##
##plt.plot(r.bins[1:], r.rdf[1:])