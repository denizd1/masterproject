from SimPEG import (
    EM, Mesh, Maps, Utils, DataMisfit, Regularization,
    Optimization, Inversion, InvProblem, Directives
)
import numpy as np
from SimPEG.EM import  FDEM, TDEM, mu_0
import matplotlib.pyplot as plt
import matplotlib
from scipy.constants import mu_0
import math
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver

sigma = 1.

csx, ncx, npadx = 5, 50, 25
csz, ncz, npadz = 5, 50, 25

hx = Utils.meshTensor(
        [(csx, ncx), (csx, npadx, 1.3)]
)
hz = Utils.meshTensor(
        [(csz, npadz, -1.3), (csz, ncz), (csz, npadz, 1.3)]
)

# define the cylindrical mesh
mesh = Mesh.CylMesh([hx, 1, hz], [0., 0., -hz.sum()/2])

rxlocs = np.array([[10, 0., 0.]])
rxtimes = np.logspace(-9, -2, 20)

#reciever 
rx = TDEM.Rx.Point_e(
locs=rxlocs, times=rxtimes, orientation='x'
)                                                                   
Aloc = np.r_[0., 0., 0.]
Bloc = np.r_[1., 0., 0.]
srcloc = np.vstack((Aloc, Bloc))
src = TDEM.Src.LineCurrent([rx], loc=srcloc, waveform=EM.TDEM.Src.StepOffWaveform())

survey = TDEM.Survey([src])


prb = TDEM.Problem3D_e(mesh,sigma=sigma) #define the problem
prb.Solver = Solver

prb.timeSteps = [(1e-06, 40), (5e-06, 40), (1e-05, 40), (5e-05, 40),
                     (0.0001, 40), (0.0005, 40)]
prb.pair(survey)

data=survey.dpred([])
rdis=9.5
xdis=9.5
efield=[]
sig_half=1
for time in range(len(rxtimes)):
        phi=math.sqrt(mu_0*sig_half/(4*rxtimes[time]))
        ee=(1/(4*np.pi*sig_half*rdis**3))*(( ( (4/math.sqrt(np.pi))*(phi**3)*(rdis**3) + (6/math.sqrt(np.pi)*phi*rdis) )*np.exp(-phi**2*rdis**2)+3*math.erfc(phi*rdis) )\
        *(xdis**2/(rdis**2))-(((4/math.sqrt(np.pi))*(phi**3)*(rdis**3)+(2/math.sqrt(np.pi))*phi*rdis))*np.exp(-phi**2*(rdis**2))+math.erfc(phi*rdis))
        efield.append(ee)

data=data-data[0]
#plot the graph
plt.semilogx(rx.times, np.absolute(data/100000))
plt.semilogx(rx.times, efield, label='analytical')
plt.legend(loc='upper left')
plt.xlabel('time (s)')
plt.ylabel('e')
plt.grid(True)
plt.show()