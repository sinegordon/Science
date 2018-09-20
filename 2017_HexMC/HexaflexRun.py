
# coding: utf-8

# # SYSTEM RUN

# In[84]:

import turtle
import numpy


rad_mul = 0.9 #float(sys.argv[1])
koch_flake = "FRFRF"
iterations = 2

for i in range(iterations):
    koch_flake = koch_flake.replace("F","FLFRFLF")

turtle.down()
coordsK = []
for move in koch_flake:
    if move == "F":
        turtle.forward(rad_mul*1.0 / (3 ** (iterations - 1)))
        coordsK.append(turtle.position())
    elif move == "L":
        turtle.left(60)
    elif move == "R":
        turtle.right(120)

coordsK = numpy.array(coordsK)
maxx = numpy.max(coordsK[:, 0])
minx = numpy.min(coordsK[:, 0])
maxy = numpy.max(coordsK[:, 1])
miny = numpy.min(coordsK[:, 1])
dx = (maxx + minx)/2.0
dy = (maxy + miny)/2.0
coordsK[:, 0] -= dx 
coordsK[:, 1] -= dy
coordsK = np.flip(coordsK, axis=0)


# In[86]:

import hoomd
import hoomd.hpmc
import ex_render
import math
import sys
import random
import numpy as np
from shapely.geometry import Polygon
hoomd.context.initialize("--mode=cpu")

rad_hex = 2.0*np.cos(np.pi/6)
vertex = []
types = []
diameters=[]
coords = []

#Hexagon1
p_hex1 = [0, 0.5*rad_hex, 0]
coords.append(p_hex1)
types.append('A')
diameters.append(2*rad_hex)
#Hexagon2
p_hex2 = [rad_hex*np.sqrt(3)/2, -rad_hex ,0]
coords.append(p_hex2)
types.append('A')
diameters.append(2*rad_hex)

print rad_mul**2

uc = hoomd.lattice.unitcell(N=2,
                            a1=[rad_hex*np.sqrt(3), 0, 0],
                            a2=[0, 3*rad_hex, 0],
                            a3=[0,   0,   1],
                            dimensions=2,
                            position=coords,
                            diameter=diameters,
                            type_name=types);

system = hoomd.init.create_lattice(unitcell=uc, n=[10, 10])
mc = hoomd.hpmc.integrate.simple_polygon(d=0.01, a=0.01, seed=42)
square_verts = coordsK
mc.shape_param.set('A', vertices=square_verts)


# In[88]:

dump_name = "koch_"+str(rad_mul**2)+".gsd"
d = hoomd.dump.gsd(dump_name, period=10, group=hoomd.group.all(), overwrite=True);


# In[89]:

hoomd.run(1001)

