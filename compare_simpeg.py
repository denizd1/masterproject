from SimPEG import (
    EM, Mesh, Maps, Utils, DataMisfit, Regularization,
    Optimization, Inversion, InvProblem, Directives
)
import numpy as np
from SimPEG.EM import  FDEM, TDEM, mu_0
import matplotlib.pyplot as plt
import matplotlib
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver


    sig_half=1  #halfspace sigma

    # Set up cylindrically symmeric mesh. 
    cs, ncx, ncz, npad = 5., 30, 10, 15
    hx = [(cs, ncx), (cs, npad, 1.3)]
    hz = [(cs, npad, -1.3), (cs, ncz), (cs, npad, 1.3)]
    mesh = Mesh.CylMesh([hx, 1, hz], '00C')
    #mesh = Mesh.CylMesh([7, 6, 10])
    mesh.plotGrid()
    plt.show()


    active = mesh.vectorCCz < 0.  #active area is under the surface
    actMap = Maps.InjectActiveCells(mesh, active, np.log(1e-8), nC=mesh.nCz)
    mapping = Maps.ExpMap(mesh) * Maps.SurjectVertical1D(mesh) * actMap  
    

    #reciever 
    rx = TDEM.Rx.Point_b(
        np.array([[50., 0., 0.]]), np.logspace(-5, -2, 20), 'z'
    )
    
    
    #define source, parameters: reciever, waveform, location
    src = TDEM.Src.MagDipole(
            [rx], waveform=TDEM.Src.StepOffWaveform(),
            loc=np.array([0., 0., 0.])
    )

    survey = TDEM.Survey([src])

    prb = TDEM.Problem3D_b(mesh, sigmaMap=mapping) #define the problem
    prb.Solver = Solver

    prb.timeSteps = [(1e-06, 40), (5e-06, 40), (1e-05, 40), (5e-05, 40),
                     (0.0001, 40), (0.0005, 40)]


    bz_ana = mu_0*EM.Analytics.hzAnalyticDipoleT(rx.locs[0][0]+1e-3,rx.times, sig_half) #analytical

    sigma = np.ones(mesh.nCz)*1e-8
    sigma[active] = sig_half
    sigma = np.log(sigma[active])
    prb.pair(survey)
    data=survey.dpred(sigma)
    fields = prb.fields(sigma)
    print(fields[src])

    #plot the graph
    plt.loglog(rx.times, abs(data), label='predicted')
    plt.loglog(rx.times, abs(bz_ana), label='analytical')
    plt.legend(loc='upper left')
    plt.xlabel('time (s)')
    plt.ylabel('bz')
    plt.grid(True)
    plt.show()

    