import matplotlib.lines as lines
from math import *
import hoomd
import hoomd.hpmc
import ex_render
import math
import sys
import random
import numpy as np

hex_rad = 1
rad_mul = 0.8
sin30 = 0.8*hex_rad*sin(pi/6)
sin60 = 0.8*hex_rad*sin(pi/3)
hex_line = [[0, hex_rad], [-sin60, sin30], [-sin60, -sin30], [0, -hex_rad], [sin60, -sin30], [sin60, sin30], 
            [0, hex_rad]]

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
        alpha = pi/3.0
        p = [(p2[0]-p1[0])/3.0, (p2[1]-p1[1])/3.0]
        p5 = [p[0]*cos(alpha) - p[1]*sin(alpha), p[1]*cos(alpha) + p[0]*sin(alpha)]
        p5 = [p5[0] + p3[0], p5[1] + p3[1]]
        new_line.append(p3)
        new_line.append(p5)
        new_line.append(p4)
    new_line.append(old_line[-1])
    new_line = np.array(new_line)
    return new_line

for i in range(2):
    hex_line = hex_flake_step(hex_line)

hoomd.context.initialize("--mode=cpu")

rad_hex = 2.0*np.cos(np.pi/6)
vertex = []
types = []
diameters=[]
coords = []
coordsK = hex_line[:-1]

#Hexagon1
p_hex1 = [0, 0.5*hex_rad, 0]
coords.append(p_hex1)
types.append('A')
diameters.append(2*hex_rad)
#Hexagon2
p_hex2 = [hex_rad*np.sqrt(3)/2, -hex_rad ,0]
coords.append(p_hex2)
types.append('A')
diameters.append(2*hex_rad)

uc = hoomd.lattice.unitcell(N=2,
                            a1=[hex_rad*np.sqrt(3), 0, 0],
                            a2=[0, 3*hex_rad, 0],
                            a3=[0,   0,   1],
                            dimensions=2,
                            position=coords,
                            diameter=diameters,
                            type_name=types);

system = hoomd.init.create_lattice(unitcell=uc, n=[10, 10])
mc = hoomd.hpmc.integrate.simple_polygon(d=0.01, a=0.01, seed=42)
square_verts = coordsK
mc.shape_param.set('A', vertices=square_verts)