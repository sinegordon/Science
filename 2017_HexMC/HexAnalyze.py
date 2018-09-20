import hoomd
import hoomd.hpmc
import ex_render
import math
import random
from scipy.spatial import distance_matrix
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from shapely.geometry import Polygon
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


trs = glob.glob("./*.gsd")
print trs
# In[6]:
k = 2
tr = ex_render.get_trajectory(trs[k])
len(tr)
vertex = []
rho = float(re.findall("\d+\.\d+", trs[k])[0])
rad_hex = np.sqrt(rho)
print rad_hex
#Hexagon1
p_hex1 = [0, 0.5*rad_hex, 0]
vertex.append(polygon(6, rad_hex, 0, p_hex1[0:2]))
#Hexagon2
p_hex2 = [rad_hex*np.sqrt(3)/2, -rad_hex ,0]
vertex.append(polygon(6, rad_hex, 0, p_hex2[0:2]))
ex_render.render_polygon_frame(tr[-1], vertex[0])

# In[ ]:

p_index = 0
box = tr[0].configuration.box
tid = tr[0].particles.typeid == p_index
N = len(tr[0].particles.position[tid, 0])
p = np.zeros((N, 3), dtype=np.float32)
sigma = 0.001
n_sigma = 10
delta = n_sigma*sigma
x_dist = [0*sigma, n_sigma*sigma, 2*n_sigma*sigma, 4*n_sigma*sigma]
x_draw = [[], [], [], []]
y_draw = [[], [], [], []]
M = 1000
frame_begin = len(tr)-1
frame_end = len(tr)
frame_period = 1
for isg in xrange(4):
    for ifr in range(frame_begin, frame_end):
        if ifr % frame_period == 0:
            print "Sigma #", isg, ": Frame #", ifr, " process"
        p[:, :] = tr[ifr].particles.position[tid, 0:3]
        p[:, 2] = 0
        p[:, 0] /= box[0]
        p[:, 1] /= box[1]
        for k in xrange(M*N):
            i = np.random.randint(0,N)
            j = np.random.randint(0,N)
            yy=(p[j,1]-p[i,1])%1.0
            if yy > 0.5: yy-=1.0
            if abs(yy) < delta:
                xx=(p[j,0]-p[i,0])%1.0
                if xx > 0.5: xx-=1.0
                if abs(xx-x_dist[isg]) < delta and i!=j:
                    x_draw[isg].append(xx/sigma)
                    y_draw[isg].append(yy/sigma)

plt.figure(figsize=(18,3.2))
plt.suptitle('Positional correlations')
plt.subplots_adjust(wspace=.3,bottom=0.15,top=0.85)
for k in xrange(4):
    plt.subplot(1,4,k+1)
    ymin=-n_sigma
    ymax=n_sigma
    xmin=x_dist[k]/sigma - n_sigma
    xmax=x_dist[k]/sigma + n_sigma
    plt.hexbin(x_draw[k],y_draw[k], gridsize=100)
    plt.axis([xmin,xmax,ymin,ymax])
    plt.xlabel("$\Delta x/\sigma$")
    plt.ylabel("$\Delta y/\sigma$",x=30.)
    plt.yticks([-n_sigma,0,n_sigma])
    plt.xticks([xmin,xmin+n_sigma,xmax])
    cb = plt.colorbar()
plt.show()


# # COMPUTE HEXATIC ORDER PARAMETER

# In[22]:


from bokeh.io import output_notebook
output_notebook()
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot
import numpy as np
import PIL.Image
import io
import IPython.display
from freud import parallel
parallel.setNumThreads(2)

# helper functions used in the notebook are below; you are free to disregard

def showarray(a, fmt='png'):
    """
    uses PIL to display an image rendered externally.
    
    Currently not used
    """
    f = io.BytesIO()
    PIL.Image.fromarray(a, mode='RGBA').save(f, fmt)
#     PIL.Image.fromarray(a, mode='RGBA').save("out.png")
    return IPython.display.display(IPython.display.Image(data=f.getvalue(), width=600))

def default_bokeh(p):
    """
    wrapper which takes the default bokeh outputs and changes them to more sensible values
    """
    p.title.text_font_size = "18pt"
    p.title.align = "center"

    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"

    p.xaxis.major_tick_in = 10
    p.xaxis.major_tick_out = 0
    p.xaxis.minor_tick_in = 5
    p.xaxis.minor_tick_out = 0

    p.yaxis.major_tick_in = 10
    p.yaxis.major_tick_out = 0
    p.yaxis.minor_tick_in = 5
    p.yaxis.minor_tick_out = 0

    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"

def cubeellipse(theta, lam=0.5, gamma=1., s=4.0, r=1., h=1.):
    """Create an RGB colormap from an input angle theta. Takes lam (a list of
    intensity values, from 0 to 1), gamma (a nonlinear weighting power),
    s (starting angle), r (number of revolutions around the circle), and
    h (a hue factor)."""
    import numpy
    lam = lam**gamma

    a = h*lam*(1 - lam)*.5
    v = numpy.array([[-.14861, 1.78277], [-.29227, -.90649], [1.97294, 0.]], dtype=numpy.float32)
    ctarray = numpy.array([numpy.cos(theta*r + s), numpy.sin(theta*r + s)], dtype=numpy.float32)
    # convert to 255 rgb
    ctarray = (lam + a*v.dot(ctarray)).T
    ctarray *= 255
    ctarray = ctarray.astype(dtype=np.int32)
    return ctarray

def local_to_global(verts, positions, orientations):
    """
    Take a list of vertices, positions, and orientations and create
    a list of vertices in the "global coordinate system" for plotting
    in bokeh
    """
    num_particles = len(positions)
    num_verts = len(verts)
    # create list of vertices in the "local reference frame" i.e.
    # centered at (0,0)
    l_verts = np.zeros(shape=(num_particles, num_verts, 2), dtype=np.float32)
    l_verts[:] = verts
    # create array of rotation matrices
    rot_mat = np.zeros(shape=(num_particles, 2, 2), dtype=np.float32)
    for i, theta in enumerate(orientations):
        rot_mat[i] = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    # rotate; uses einsum for speed; please see numpy documentation
    # for more information
    r_verts = np.einsum("lij,lkj->lki", rot_mat, l_verts)
    # now translate to global coordinates
    # need to create a position array with same shape as vertex array
    l_pos = np.zeros(shape=(num_particles, num_verts, 2), dtype=np.float32)
    for i in range(num_particles):
        for j in range(len(verts)):
            l_pos[i,j] = positions[i]
    # translate
    output_array = np.add(r_verts, l_pos)
    return output_array

def clamp(x):
    """
    limit values between 0 and 255
    http://stackoverflow.com/questions/3380726/converting-a-rgb-color-tuple-to-a-six-digit-code-in-python
    """
    return max(0, min(x, 255))


# In[23]:


from freud import box, order

def GetHexatic(frame_ind, rmax):
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


# render in bokeh
# vertex positions for hexagons
data = GetHexatic(-1,1.2)
psi_k = data[0]
avg_psi_k = data[1]
ang_data = data[2]
pos_data = data[3]
verts = vertex[0]
verts = np.array(verts)
# create array of transformed positions
patches = local_to_global(verts, pos_data[:, 0:2], ang_data)
# create an array of angles relative to the average
a = np.angle(psi_k) - np.angle(avg_psi_k)
# turn into an rgb array of tuples
color = [tuple(cubeellipse(x)) for x in a]
# bokeh (as of this version) requires hex colors, so convert rgb to hex
hex_color = ["#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b)) for (r,g,b) in color]
# plot
p = figure(title="Hexatic Order Parameter visualization")
p.patches(xs=patches[:,:,0].tolist(), ys=patches[:,:,1].tolist(),
    fill_color=hex_color, line_color="black")
default_bokeh(p)
show(p)

# # COMPUTE  $G_6(r)$

# In[112]:

data = GetHexatic(-1, 4)
psi_k = data[0]
avg_psi_k = data[1]
ang_data = data[2]
pos_data = data[3]
g6ij_mas = np.outer(psi_k, np.conj(psi_k))
rij_mas = distance_matrix(pos_data[:, 0:2], pos_data[:, 0:2])
dr = 0.5
r_mas = []
g6_mas = []
for r in np.arange(1.0, 10, 0.1):
    inds = (rij_mas > r - dr)*(rij_mas < r + dr)*(rij_mas > 0.0)
    g6_mas.append(np.mean(g6ij_mas[inds]))
    r_mas.append(r)
g6_mas = np.abs(np.array(g6_mas))
r_mas = np.array(r_mas)


# In[124]:

plt.plot(r_mas, np.abs(g6_mas))
plt.show()


# In[ ]:




