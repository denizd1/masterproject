from SimPEG import (
    EM, Mesh, Maps, Utils, DataMisfit, Regularization,
    Optimization, Inversion, InvProblem, Directives
)
import numpy as np
from SimPEG.EM import  FDEM, TDEM, mu_0
import matplotlib.pyplot as plt
import matplotlib
from scipy.special import erfc, erf
from scipy.constants import pi, epsilon_0
import math
from pymatsolver import Pardiso as Solver

sig_half=1.
# Cell sizes
cs, nc, npad = 10., 10, 10
hx = [(cs, npad, -1.3), (2, 25), (cs, npad, 1.3)]
hy = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
hz = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
mesh = Mesh.TensorMesh([hx, hy, hz], 'CCC')
# nc, npad = 5, 5
# cs = 5
# dfine = 100
# nc = dfine/cs

# hx = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
# hy = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
# hz = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
# mesh = Mesh.TensorMesh([hx, hy, hz], 'CCC')
mesh.plotGrid()
plt.show()
print(mesh)

rxlocs = np.array([[30, 0., 0.]])
rxtimes = np.logspace(-5, -3, 31)
#reciever 
rx = TDEM.Rx.Point_e(
    locs=rxlocs, times=rxtimes, orientation='x'
)
Aloc = np.r_[-0.5, 0., 0.]
Bloc = np.r_[0.5, 0., 0.]
srcloc = np.vstack((Aloc, Bloc))


#define source, parameters: reciever, waveform, location
src = TDEM.Src.LineCurrent([rx], loc=srcloc, waveform=EM.TDEM.Src.StepOffWaveform())

survey = TDEM.Survey([src])

prb = TDEM.Problem3D_e(mesh, sigma=sig_half) 
prb.Solver = Solver

prb.timeSteps = [(1e-06, 20), (1e-05, 20), (0.0001, 20)]


prb.pair(survey)
data=survey.dpred([])

rdis=np.sqrt( 30**2. + 0**2. + 0**2.)
xdis=30


theta=np.sqrt(mu_0*sig_half/(4*rxtimes))
errfunc=erfc(theta*rdis)
ee=(1/(4*np.pi*sig_half*rdis**3))*(( ( (4/np.sqrt(np.pi))*(theta**3)*(rdis**3) + (6/np.sqrt(np.pi)*theta*rdis) )*np.exp(-theta**2*rdis**2)+3*errfunc )\
*(xdis**2/(rdis**2))-((((4/np.sqrt(np.pi))*(theta**3)*(rdis**3)+(2/np.sqrt(np.pi))*theta*rdis))*np.exp(-theta**2*(rdis**2))+errfunc))


data=data-data[0]
plt.semilogx(rx.times, -data)
plt.semilogx(rx.times, ee, label='analytical')
plt.legend(loc='upper left')
plt.xlabel('time (s)')
plt.ylabel('e')
plt.grid(True)
plt.show()
