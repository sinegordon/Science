import hoomd
import hoomd.hpmc
import ex_render
import math
import random
from scipy.spatial import distance_matrix
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from freud import parallel, box, density, order
import pickle
parallel.setNumThreads(2)
hoomd.context.initialize("--mode=cpu")

def polygon(n, a, alpha0, p0):
    l = []
    p = [0, a]
    p = [p[0]*math.cos(alpha0) - p[1]*math.sin(alpha0), p[1]*math.cos(alpha0) + p[0]*math.sin(alpha0)]
    alpha = 2*math.pi/n
    for i in range(n):
        p = [p[0]*math.cos(alpha) - p[1]*math.sin(alpha), p[1]*math.cos(alpha) + p[0]*math.sin(alpha)]
        l.append([p[0]+p0[0], p[1]+p0[1]])
    return l

def hex_flake_step(old_line):
    new_line = []
    for i in range(len(old_line)-1):
        p1 = old_line[i]
        p2 = old_line[i+1]
        dpx = (p2[0] - p1[0])/3.0
        dpy = (p2[1] - p1[1])/3.0
        new_line.append(p1)
        p3 = [p1[0]+dpx, p1[1]+dpy]
        p4 = [p1[0]+2*dpx, p1[1]+2*dpy]
        alpha = math.pi/3.0
        p = [(p2[0]-p1[0])/3.0, (p2[1]-p1[1])/3.0]
        p5 = [p[0]*math.cos(alpha) - p[1]*math.sin(alpha), p[1]*math.cos(alpha) + p[0]*math.sin(alpha)]
        p5 = [p5[0] + p3[0], p5[1] + p3[1]]
        new_line.append(p3)
        new_line.append(p5)
        new_line.append(p4)
    p1 = old_line[-1]
    p2 = old_line[0]
    dpx = (p2[0] - p1[0])/3.0
    dpy = (p2[1] - p1[1])/3.0
    new_line.append(p1)
    p3 = [p1[0]+dpx, p1[1]+dpy]
    p4 = [p1[0]+2*dpx, p1[1]+2*dpy]
    alpha = math.pi/3.0
    p = [(p2[0]-p1[0])/3.0, (p2[1]-p1[1])/3.0]
    p5 = [p[0]*math.cos(alpha) - p[1]*math.sin(alpha), p[1]*math.cos(alpha) + p[0]*math.sin(alpha)]
    p5 = [p5[0] + p3[0], p5[1] + p3[1]]
    new_line.append(p3)
    new_line.append(p5)
    new_line.append(p4)
    #new_line.append(old_line[-1])
    new_line = np.array(new_line)
    return new_line

def GetHexatic(tr, frame_ind, rmax):
    # create hexatic object
    hex_order = order.HexOrderParameter(rmax=rmax, k=6, n=6);
    p_index = 0
    tid = tr[frame_ind].particles.typeid == p_index
    n_par = len(tr[frame_ind].particles.position[tid, 0])
    pos_data = np.zeros((n_par, 3), dtype=np.float32)
    # compute the hexatic order
    pos_data[:, :] = tr[frame_ind].particles.position[tid, 0:3]
    pos_data[:, 2] = 0
    ang_data = 2*np.arctan2(np.copy(tr[frame_ind].particles.orientation[tid,3]), 
                             np.copy(tr[frame_ind].particles.orientation[tid,0]))
    # create box
    fbox = box.Box(Lx=tr[0].configuration.box[0], Ly=tr[0].configuration.box[1], is2D=True)
    # compute hexatic order for 6 nearest neighbors
    hex_order.compute(fbox, pos_data)
    # get values from freud object
    psi_k = hex_order.getPsi()
    avg_psi_k = np.mean(psi_k)
    return (psi_k, avg_psi_k, ang_data, pos_data)

def GetAll(fname):
    #RETURN [HEX_AREA, HEX_CELL_AREA, PHI, RMAS1, RDF, RMAS2, G6MAS, MSD]
    ret = []
    
    #LOAD TRAJECTORY
    tr = ex_render.get_trajectory(fname)
    print "LOAD " + fname + " DONE"
    vertex = []
    rad_hex = float(re.findall("\d+\.\d+", fname)[0])
    rho = rad_hex**2
    #Hexagon1
    p_hex1 = [0, 0.5*rad_hex, 0]
    hex1_line = polygon(6, rad_hex, 0, p_hex1[0:2])
    hex_line = np.copy(hex1_line)
    it = int(re.findall("\d", fname)[-1])
    for i in range(it):
        hex_line = hex_flake_step(hex_line)
    vertex.append(hex1_line)
    hex1 = Polygon(hex_line)
    #Cell hexagon
    hex_cell_line = polygon(6, 1, 0, p_hex1[0:2])
    hex_cell = Polygon(hex_cell_line)
    ret.append(hex1.area)
    ret.append(hex_cell.area)
    ret.append(hex1.area/hex_cell.area)
    #Hexagon2
    p_hex2 = [rad_hex*np.sqrt(3)/2, -rad_hex ,0]
    vertex.append(polygon(6, rad_hex, 0, p_hex2[0:2]))

    #COMPUTE 1RDF
    print "COMPUTE RDF"
    p_index = 0
    tid = tr[0].particles.typeid == p_index
    n_par = len(tr[0].particles.position[tid, 0])
    p = np.zeros((n_par, 3), dtype=np.float32)
    rdf = density.RDF(rmax=10, dr=0.1)
    # compute the rdf for for all frames except the first (your syntax will vary based on your reader)
    frame_begin = 1
    frame_end = len(tr)
    frame_period = 1
    # compute the rdf for for first frame
    fbox = box.Box(Lx=tr[0].configuration.box[0], Ly=tr[0].configuration.box[1], is2D=True)
    rdf.compute(fbox, p, p)
    p[:, :] = tr[frame_begin].particles.position[tid, 0:3]
    p[:, 2] = 0
    for i in range(frame_begin, frame_end):
        # read box, position data
        #if i % frame_period == 0:
        #    print "Frame #", i, " process"
        p[:, :] = tr[i].particles.position[tid, 0:3]
        p[:, 2] = 0
        # create the freud box object
        fbox = box.Box(Lx=tr[0].configuration.box[0], Ly=tr[0].configuration.box[1], is2D=True)
        # accumulate
        rdf.accumulate(fbox, p, p)
    ret.append(rdf.getR())
    ret.append(rdf.getRDF())
    
    # COMPUTE G6
    print "COMPUTE G6"
    data = GetHexatic(tr, -1, 5)
    k = 4
    psi_k = data[0][0::k]
    avg_psi_k = data[1]
    ang_data = data[2]
    pos_data = data[3][0::k]
    g6ij_mas = np.outer(psi_k, np.conj(psi_k))
    rij_mas = distance_matrix(pos_data[:, 0:2], pos_data[:, 0:2])
    dr = 0.2
    r_mas = []
    g6_mas = []
    for r in np.arange(1.0, 20, 0.1):
        inds = (rij_mas > r - dr)*(rij_mas < r + dr)*(rij_mas > 0.0)
        g6_mas.append(np.mean(g6ij_mas[inds]))
        r_mas.append(r)
    ret.append(r_mas)
    ret.append(g6_mas)
    
    #COMPUTE MSD
    print "COMPUTE MSD"
    p_index = 0
    tid = tr[0].particles.typeid == p_index
    n_par = len(tr[0].particles.position[tid, 0])
    pos0 = tr[0].particles.position[tid, 0:2]
    Lx=tr[0].configuration.box[0]
    Ly=tr[0].configuration.box[1]
    frame_begin = 1
    frame_end = len(tr)
    frame_period = 1
    msd = []
    for i in range(frame_begin, frame_end):
        pos = tr[i].particles.position[tid, 0:2]
        dx = pos[:, 0] - pos0[:, 0]
        dy = pos[:, 1] - pos0[:, 1]
        indsx = np.abs(dx) < Lx/2
        indsy = np.abs(dy) < Ly/2
        inds = indsx * indsy
        msd.append(np.mean(dy[inds]**2 + dx[inds]**2))
    ret.append(msd)
    return ret

trs = glob.glob("/Users/user/Documents/Work/2017_HexMC/*.gsd")
len(trs)
print trs

mas = {}
for fname in trs:
    out = GetAll(fname)
    mas[fname] = out

f = open("./data.pickle","wb")
pickle.dump(mas, f)
f.close()