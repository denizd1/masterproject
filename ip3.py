from __future__ import print_function
import unittest
from SimPEG import Mesh, Utils, Maps
import numpy as np
import SimPEG.EM.Static.DC as DC
import SimPEG.EM.Static.IP as IP
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver
import matplotlib.pyplot as plt


cs = 12.5
hx = [(cs, 7, -1.3), (cs, 61), (cs, 7, 1.3)]
hy = [(cs, 7, -1.3), (cs, 20)]
mesh = Mesh.TensorMesh([hx, hy], x0="CN")

x = np.linspace(-200, 200., 20)
M = Utils.ndgrid(x-12.5, np.r_[0.])
N = Utils.ndgrid(x+12.5, np.r_[0.])
A0loc = np.r_[-150, 0.]
A1loc = np.r_[-130, 0.]
B0loc = np.r_[-130, 0.]
B1loc = np.r_[-110, 0.]

rx = DC.Rx.Dipole_ky(M, N)
src0 = DC.Src.Pole([rx], A0loc)

src0_ip = DC.Src.Pole([rx], A0loc)


srcLists_ip = [src0_ip]
surveyDC = DC.Survey_ky([src0])

sigmaInf = np.ones(mesh.nC) * 1.
blkind = Utils.ModelBuilder.getIndicesSphere(
    np.r_[0, -150], 40, mesh.gridCC)

eta = np.zeros(mesh.nC)
eta[blkind] = 0.1
sigma0 = sigmaInf * (1.-eta)
# mesh.plotImage(sigma0)
# plt.show()

problemDC = DC.Problem2D_N(mesh, sigmaMap=Maps.IdentityMap(mesh))
problemDC.Solver = Solver
problemDC.pair(surveyDC)
data0 = surveyDC.dpred(sigma0)
datainf = surveyDC.dpred(sigmaInf)
problemIP = IP.Problem2D_N(
    mesh,
    sigma=sigmaInf,
    etaMap=Maps.IdentityMap(mesh),
)
problemIP.Solver = Solver
surveyIP = IP.Survey(srcLists_ip)
problemIP.pair(surveyIP)
data_full = data0 - datainf
data = surveyIP.dpred(eta)
plt.plot(data_full)
plt.plot(data, 'k.')
plt.show()