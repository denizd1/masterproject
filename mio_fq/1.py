import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from SimPEG import Mesh, EM, Utils
from simpegEMIP.FDEM import Problem3D_e
from SimPEG.EM import FDEM
from SimPEG.EM.FDEM.SrcFDEM import LineCurrent
from SimPEG import Mesh
import numpy as np
from pymatsolver import Pardiso as Solver
import time
import os
from tqdm import tqdm

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
datafile=os.path.join(THIS_FOLDER,'freqs.txt')

core_domain_x = np.r_[-250., 250.]  # extent of uniform cells in the x-direction
core_domain_z = np.r_[-200., 0.]  # extent of uniform cells in the z-direction
core_domain_y = np.r_[-250.,250.]
csx, csy, csz = 10, 15, 15

# Number of core cells in each direction
ncy = int(np.diff(core_domain_y)/csy)
ncx = int(np.diff(core_domain_x)/csx)
ncz = int(np.diff(core_domain_z)/csz)
# Number of padding cells to add in each direction
npad = 10
# Vectors of cell lengths in each direction with padding
hx = [(csx, npad, -1.5), (csx, ncx), (csx, npad, 1.5)]
hy = [(csy, npad, -1.5), (csy, ncy), (csy, npad, 1.5)]
hz = [(csz, npad, -1.5), (csz, ncz), (csz, npad, 1.5)]
# Create mesh and center it
mesh = Mesh.TensorMesh([hx, hy, hz])
#print(-mesh.hx.sum()/2)
mesh.x0 = np.r_[
    -mesh.hx.sum()/2, -mesh.hy.sum()/2., -mesh.hz[:npad+ncz].sum()
]


p0 = np.array([-85, -85, -50])
p1 = np.array([85, 85, -100])
inds = Utils.ModelBuilder.getIndicesBlock(p0, p1, mesh.gridCC)
sigmaInf = np.ones(mesh.nC) *1e-8
sigmaInf[mesh.gridCC[:,2]<0.] = 1
eta = np.ones(mesh.nC) * 0
tau = np.ones(mesh.nC) *1e-3
c = np.ones(mesh.nC)

#eta[mesh.gridCC[:,2]<0.] = 0.1
#tau[mesh.gridCC[:,2]<0.] = 0.1
#c[mesh.gridCC[:,2]<0.] = 0.5
sigmaInf[inds]=0.1
eta[inds] = 0.1
tau[inds] = 1
c[inds] = 0.5

# out = mesh.plotSlice(sigmaInf, grid=True, normal='X')
# plt.colorbar(out[0])
# plt.title("Conductivity (S/m)")
# plt.gca().set_aspect(1)
# plt.show()

# out = mesh.plotSlice(eta, grid=True, normal='X')
# plt.colorbar(out[0])
# plt.title("Chargeability (V/V)")
# plt.gca().set_aspect(1)
# plt.show()

x = np.linspace(-200, 200, 21)
y = np.linspace(-200, 200, 21)
ix = np.arange(len(x))
iy = np.arange(len(y))
rx_locations = Utils.ndgrid(x, y, np.r_[0.])
src_locations = np.array(
    [[-0.5, 0, 0],[0.5, 0, 0]]
)



# plt.plot(rx_locations[:,0], rx_locations[:,1], 'k.')
# plt.plot(src_locations[:,0], src_locations[:,1], 'r-')
# ir = np.arange(len(rx_locations))
# for (_x, _y, _z), i in zip(rx_locations, ir):
#     plt.text(_x, _y, str(i),horizontalalignment='center')
# # plt.xticks(ticks=x,labels=ix)
# # plt.yticks(ticks=y,labels=iy)
# plt.gca().set_aspect(1)
# plt.show()


rx_r = FDEM.Rx.Point_e(rx_locations, orientation='x', component='real')
rx_i = FDEM.Rx.Point_e(rx_locations, orientation='x', component='imag')

with open(datafile, "r") as f:
     f = np.array([float(line) for line in f])


freqs=f[:20]/(2*np.pi)
#freqs=np.logspace(-3,3,100)
dataarray=[]
for f in tqdm(range(len(freqs))):
    src = LineCurrent([rx_r, rx_i], freq=freqs[f], loc=src_locations)
    survey = FDEM.Survey([src])
    problem = Problem3D_e(
        mesh, 
        sigmaInf=sigmaInf, eta=eta, tau=tau, c=c,
        Solver=Solver
    )
    problem.pair(survey)

    data = survey.dpred([])
    d=data.reshape((441, 2),order='F')[210:231,:]
    dataarray.append(d.tolist())

#print(time.time()-a)
#print(dataarray)
np.save('t1_f1_hs.npy',dataarray)
