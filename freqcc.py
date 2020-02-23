from SimPEG import (
    EM, Mesh, Maps, Utils, DataMisfit, Regularization,
    Optimization, Inversion, InvProblem, Directives
)
import numpy as np
from SimPEG.EM import  FDEM, TDEM, mu_0
import matplotlib.pyplot as plt
import matplotlib
from scipy.constants import pi, epsilon_0
import math
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver

    sig_half=1 +2j
    cs, nc, npad = 5., 5, 5
    hx = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
    hy = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
    hz = [(cs, npad, -1.3), (cs, nc), (cs, npad, 1.3)]
    mesh = Mesh.TensorMesh([hx, hy, hz], 'CCC')
    sigmaInf = np.ones(mesh.nC) * 0.001+1j
    layerind = (np.logical_and(mesh.gridCC[:,2]<0, mesh.gridCC[:,2]>-50.)) & (np.logical_and(mesh.gridCC[:,0]<500, mesh.gridCC[:,0]>100.))

    active = mesh.vectorCCz < 0.
    actMap = Maps.InjectActiveCells(mesh, active, np.log(1e8), nC=mesh.nCz)
    mapping = Maps.ExpMap(mesh) * Maps.SurjectVertical1D(mesh) * actMap
    sigmaInf[~layerind] = sig_half


    
    rxlocs = np.array([[50, 0., 0.]])
    rxtimes = np.logspace(-9, -2, 20)
    
    #reciever 
    rx = FDEM.Rx.Point_e(
        locs=rxlocs, orientation='x', component='real'
    )
    Aloc = np.r_[0., 0., 0.]
    Bloc = np.r_[1, 0., 0.]
    srcloc = np.vstack((Aloc, Bloc))
    
    
    #define source, parameters: reciever, waveform, location
    src = FDEM.Src.RawVec_e([rx], freq=1, s_e=srcloc)

    survey = FDEM.Survey([src])

    prb = FDEM.Problem3D_e(mesh, sigmaMap=mapping) 
    prb.Solver = Solver

    prb.timeSteps = [(1e-06, 40), (5e-06, 40), (1e-05, 40), (5e-05, 40),
                     (0.0001, 40), (0.0005, 40)]




    prb.pair(survey)
    data=survey.dpred(sigmaInf)


    