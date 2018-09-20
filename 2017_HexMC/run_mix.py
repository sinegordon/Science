from math import *
import hoomd
import hoomd.hpmc
import numpy as np
import sys

hex_rad = 1.0
rad_mul = float(sys.argv[1])
iters_count1 = int(sys.argv[2])
iters_count2 = int(sys.argv[3])
if len(sys.argv) == 5:
    gpuid = sys.argv[4]
sin30 = rad_mul*hex_rad*sin(pi/6)
sin60 = rad_mul*hex_rad*sin(pi/3)
hex_line1 = [[0, rad_mul*hex_rad], [-sin60, sin30], [-sin60, -sin30], [0, -rad_mul*hex_rad], [sin60, -sin30], [sin60, sin30], [0, rad_mul*hex_rad]]
hex_line2 = [[0, rad_mul*hex_rad], [-sin60, sin30], [-sin60, -sin30], [0, -rad_mul*hex_rad], [sin60, -sin30], [sin60, sin30], [0, rad_mul*hex_rad]]

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
    p1 = old_line[-1]
    p2 = old_line[0]
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
    #new_line.append(old_line[-1])
    new_line = np.array(new_line)
    return new_line

for i in range(iters_count1):
    hex_line1 = hex_flake_step(hex_line1)
for i in range(iters_count2):
    hex_line2 = hex_flake_step(hex_line2)
if len(sys.argv) == 5:
    hoomd.context.initialize("--mode=gpu --gpu="+gpuid)
else:
    hoomd.context.initialize("--mode=gpu")

rad_hex = 2.0*np.cos(np.pi/6)
vertex = []
types = []
diameters=[]
coords = []
coordsK1 = hex_line1[:-1]
coordsK2 = hex_line2[:-1]

#Hexagon1
p_hex1 = [0, 0.5*hex_rad, 0]
coords.append(p_hex1)
types.append('A')
diameters.append(2*hex_rad)
#Hexagon2
p_hex2 = [hex_rad*np.sqrt(3)/2, -hex_rad ,0]
coords.append(p_hex2)
types.append('B')
diameters.append(2*hex_rad)
uc = hoomd.lattice.unitcell(N=2,
                            a1=[hex_rad*np.sqrt(3), 0, 0],
                            a2=[0, 3*hex_rad, 0],
                            a3=[0,   0,   1],
                            dimensions=2,
                            position=coords,
                            diameter=diameters,
                            type_name=types);

system = hoomd.init.create_lattice(unitcell=uc, n=[128, 128])
if iters_count1 == 0 and iters_count2 == 0:
	mc = hoomd.hpmc.integrate.convex_polygon(d=0.01, a=0.01, seed=4636722)
else:
	mc = hoomd.hpmc.integrate.simple_polygon(d=0.01, a=0.01, seed=4633252)
hex_verts1 = coordsK1
hex_verts2 = coordsK2
mc.shape_param.set('A', vertices=hex_verts1)
mc.shape_param.set('B', vertices=hex_verts2)
hoomd.run(9000000)
d = hoomd.dump.gsd("hex_flake_"+str(rad_mul)+"_"+str(iters_count)+".gsd", period=10000, group=hoomd.group.all(), overwrite=True)
hoomd.run(1000000)
