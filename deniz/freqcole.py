import matplotlib.pyplot as plt
import numpy as np
from SimPEG import Mesh, EM, Utils
from simpegEMIP.FDEM import Problem3D_e
from SimPEG.EM import FDEM
from SimPEG.EM.FDEM.SrcFDEM import LineCurrent
from SimPEG import Mesh
import numpy as np
from pymatsolver import Pardiso as Solver

csx, csy, csz = 20., 20., 35
# Number of core cells in each direction
ncx, ncy, ncz = 25, 25, 20
# Number of padding cells to add in each direction
npad = 7
# Vectors of cell lengths in each direction with padding
hx = [(csx, npad, -1.5), (csx, ncx), (csx, npad, 1.5)]
hy = [(csy, npad, -1.5), (csy, ncy), (csy, npad, 1.5)]
hz = [(csz, 0, -1.5), (csz, ncz), (csy, 5, 1.5)]
mesh = Mesh.TensorMesh([hx, hy, hz], x0="CCC")
print(mesh)
mesh.plotGrid()
plt.show()

p0 = np.array([-300, -200, -50])
p1 = np.array([300, 200, -300])
inds = Utils.ModelBuilder.getIndicesBlock(p0, p1, mesh.gridCC)
sigmaInf = np.ones(mesh.nC) * 1
sigmaInf[mesh.gridCC[:,2]<0.] = 1.
eta = np.ones(mesh.nC) * 0
tau = np.ones(mesh.nC) * 0
c = np.ones(mesh.nC)
sigmaInf[inds]=1
eta[inds] = 0
tau[inds] = 0
c[inds] = 0

out = mesh.plotSlice(sigmaInf, grid=True, normal='X')
plt.colorbar(out[0])
plt.title("Conductivity (S/m)")
plt.gca().set_aspect(1)
plt.show()

out = mesh.plotSlice(eta, grid=True, normal='X')
plt.colorbar(out[0])
plt.title("Chargeability (V/V)")
plt.gca().set_aspect(1)
plt.show()

x = np.linspace(-400, 400, 21)
y = np.linspace(-400, 400, 21)
ix = np.arange(len(x))
iy = np.arange(len(y))
rx_locations = Utils.ndgrid(x, y, np.r_[0.])
src_locations = np.array(
    [[-0.5, 0, 0],[0.5, 0, 0]]
)



plt.plot(rx_locations[:,0], rx_locations[:,1], 'k.')
plt.plot(src_locations[:,0], src_locations[:,1], 'r-')
ir = np.arange(len(rx_locations))
for (_x, _y, _z), i in zip(rx_locations, ir):
    plt.text(_x, _y, str(i),horizontalalignment='center')
# plt.xticks(ticks=x,labels=ix)
# plt.yticks(ticks=y,labels=iy)
plt.gca().set_aspect(1)
plt.show()


rx_r = FDEM.Rx.Point_e(rx_locations, orientation='x', component='real')
rx_i = FDEM.Rx.Point_e(rx_locations, orientation='x', component='imag')
freqs=np.logspace(-3,3,25)

dataarray=[]
for k in range(len(freqs)):
    src = LineCurrent([rx_r, rx_i], freq=freqs[k], loc=src_locations)
    survey = FDEM.Survey([src])
    problem = Problem3D_e(
        mesh, 
        sigmaInf=sigmaInf, eta=eta, tau=tau, c=c,
        Solver=Solver
    )
    problem.pair(survey)

    data = survey.dpred([])
    d=data.reshape((441, 2),order='F')[225,:]
    dataarray.append(d.tolist())

np.save('wholespace.npy',dataarray)


#Utils.plot2Ddata(rx_locations, data.reshape((441, 2), order='F')[:,0])
#plt.show()
#Utils.plot2Ddata(rx_locations, data.reshape((441, 2), order='F')[:,1])
#plt.show()