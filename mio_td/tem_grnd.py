from SimPEG import Mesh, EM, Utils
import numpy as np
from simpegEMIP.TDEM.ProblemTDEMIP import Problem3D_e
from simpegEMIP.TDEM import Survey
from simpegEMIP.TDEM.Rx import Point_e
from pymatsolver import Pardiso




#core_domain_x = np.r_[-1000., 500.]  # extent of uniform cells in the x-direction
core_domain_x = np.r_[-500.,500.]
core_domain_z = np.r_[-400., 0.]  # extent of uniform cells in the z-direction
core_domain_y = np.r_[-450.,450.]
csx, csy, csz = 10, 15, 10

# Number of core cells in each direction
ncy = 10
ncx = 10
ncz = 10
# Number of padding cells to add in each direction
npad = 10
#Vectors of cell lengths in each direction with padding
hx = [(csx, npad, -1.5), (csx, ncx), (csx, npad, 1.5)]
hy = [(csy, npad, -1.5), (csy, ncy), (csy, npad, 1.5)]
hz = [(csz, npad, -1.5), (csz, ncz), (csz, npad, 1.5)]
# Create mesh and center it
mesh = Mesh.TensorMesh([hx, hy, hz])
mesh.x0 = np.r_[
    -mesh.hx.sum()/2, -mesh.hy.sum()/2., -mesh.hz[:npad+ncz].sum()
]

def Circle2D(xc, r, n):
    theta = np.linspace(-np.pi, np.pi, n)
    x = r*np.cos(Utils.mkvc(theta))+xc[0]
    y = r*np.sin(Utils.mkvc(theta))+xc[1]

    return np.c_[x, y]

x0 = np.r_[-450, -450, 0.] 
x1 = np.r_[450, -450, 0.]
x2 = np.r_[450, 450, 0.]
x3 = np.r_[-450, 450, 0.]
#pts0 = np.vstack((x0, x1, x2, x3, x0))
pts0=np.array([x0,x1,x2,x3,x0])
x0 = np.r_[-250, -250, -600.]
x1 = np.r_[250, -250, -600.]
x2 = np.r_[250, 250, -600.]
x3 = np.r_[-250, 250, -600.]
#pts1 = np.vstack((x0, x1, x2, x3, x0))
pts1=np.array([x0,x1,x2,x3,x0])
#pts = np.vstack((pts0, pts1))
pts = np.concatenate((pts0, pts1))
inds_porphyry = Utils.ModelBuilder.PolygonInd(mesh, pts)
inds_overburden = np.logical_and(mesh.gridCC[:,2]<0., mesh.gridCC[:,2]>-30.)
