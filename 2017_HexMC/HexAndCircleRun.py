
# coding: utf-8

import hoomd
import hoomd.hpmc
import ex_render
import math
import random
from shapely.geometry import Polygon
import sys
hoomd.context.initialize("--mode=cpu");


def polygon(n, a, alpha0, p0):
    l = []
    p = [0, a]
    p = [p[0]*math.cos(alpha0) - p[1]*math.sin(alpha0), p[1]*math.cos(alpha0) + p[0]*math.sin(alpha0)]
    alpha = 2*math.pi/n
    for i in range(n):
        p = [p[0]*math.cos(alpha) - p[1]*math.sin(alpha), p[1]*math.cos(alpha) + p[0]*math.sin(alpha)]
        l.append([p[0]+p0[0], p[1]+p0[1]])
    return l

mult_poly = float(sys.argv[1])
count_box = 6
dx, dy = 1.0/count_box, 1.0/count_box
rad_poly = mult_poly/(2*count_box)
rad_hex = 0.35
n_vertex = 21
coords = []
vertex = []
types = []
diameters=[]

#Hexagon
p_hex = [0,0,0]
coords.append(p_hex)
vertex.append(polygon(6, rad_hex, math.pi/6, p_hex[0:2]))
types.append('A')
diameters.append(2*rad_hex)
hex_coords = []
for v in vertex[0]:
    hex_coords.append((p_hex[0]+v[0], p_hex[1]+v[1]))
hex_obj = Polygon(hex_coords)

#Small-polygons (disks) lattice constructor
count_poly = 0
for i in range(count_box):
    for j in range(count_box):
        x, y = -0.5 + dx/2 + i*dx, -0.5 + dy/2 + j*dy
        poly_coords = polygon(n_vertex, rad_poly, 0, [x, y])
        poly_obj = Polygon(poly_coords)
        #Adding disk to lattice if the one dont intersecting hexagone
        if hex_obj.intersects(poly_obj) == False:
            coords.append([x, y, 0])
            vertex.append(polygon(n_vertex, rad_poly, 0, [0, 0]))
            types.append('B')
            diameters.append(2*rad_poly)
            count_poly += 1
sqr_hex = 6*rad_hex**2*math.cos(math.pi/6)*math.sin(math.pi/6)
sqr_poly = count_poly*n_vertex*rad_poly**2*math.cos(math.pi/n_vertex)*math.sin(math.pi/n_vertex)
etta = (sqr_hex + sqr_poly)
print etta


uc = hoomd.lattice.unitcell(N=count_poly + 1,
                            a1=[1.0, 0, 0],
                            a2=[0, 1.0, 0],
                            a3=[0, 0, 1.0],
                            dimensions=2,
                            position=coords,
                            diameter=diameters,
                            type_name=types);


system = hoomd.init.create_lattice(unitcell=uc, n=[100, 100]);

mc = hoomd.hpmc.integrate.convex_polygon(d=0.01, a=0.01, seed=189)
square_vertsA = vertex[0]
mc.shape_param.set('A', vertices=square_vertsA)
square_vertsB = vertex[1]
mc.shape_param.set('B', vertices=square_vertsB)

#ex_render.render_2polygons_frame(system.take_snapshot(all=True), square_vertsA, square_vertsB)

d = hoomd.dump.gsd("trajectory_" + str(etta) + ".gsd", period=10000, group=hoomd.group.all(), overwrite=True)

hoomd.run(101000)
#tr = ex_render.get_trajectory("trajectory.gsd")

#ex_render.render_2polygons_frame(tr[-1], square_vertsA, square_vertsB)