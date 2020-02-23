from SimPEG import (Mesh, Utils, Maps, DataMisfit,
                    Regularization, Optimization, Inversion, InvProblem, Tests)
import numpy as np
from SimPEG.EM.Static import SIP
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver
import matplotlib.pyplot as plt
import matplotlib

cs = 25.
hx = [(cs, 0, -1.3), (cs, 21), (cs, 0, 1.3)]
hz = [(cs, 0, -1.3), (cs, 20)]
mesh = Mesh.TensorMesh([hx, hz], x0="CN")
blkind0 = Utils.ModelBuilder.getIndicesSphere(
    np.r_[-100., -200.], 75., mesh.gridCC
)
blkind1 = Utils.ModelBuilder.getIndicesSphere(
    np.r_[100., -200.], 75., mesh.gridCC
)

sigma = np.ones(mesh.nC) * 1e-2
eta = np.zeros(mesh.nC)
tau = np.ones_like(sigma) * 1.
eta[blkind0] = 0.1
eta[blkind1] = 0.1
tau[blkind0] = 0.1
tau[blkind1] = 0.1

x = mesh.vectorCCx[(mesh.vectorCCx > -155.) & (mesh.vectorCCx < 155.)]

# mesh.plotImage(eta)
# plt.show()

Aloc = np.r_[-200., 0.]
Bloc = np.r_[200., 0.]
M = Utils.ndgrid(x-25., np.r_[0.])
N = Utils.ndgrid(x+25., np.r_[0.])

times = np.arange(10)*1e-3 + 1e-3
rx = SIP.Rx.Dipole(M, N, times)
src = SIP.Src.Pole([rx], Aloc)
survey = SIP.Survey([rx],[src])
wires = Maps.Wires(('eta', mesh.nC), ('taui', mesh.nC))
problem = SIP.Problem2D_N(
    mesh,
    sigma=sigma,
    etaMap=wires.eta,
    tauiMap=wires.taui,
    verbose = False
)
problem.Solver = Solver
problem.pair(survey)
mSynth = np.r_[eta, 1./tau]
problem.model = mSynth
data=survey.dpred(problem.model)
plt.semilogx(data)
plt.show()
