from SimPEG import *
from SimPEG import EM
from scipy.constants import mu_0
import numpy as np
import scipy.sparse as sp
from simpegEMIP.TDEM import geteref, Problem3DIP_Linear, SurveyLinear
from simpegEMIP.TDEM import Survey, Rx
from simpegEMIP.TDEM import Problem3DEM_e, Problem3D_e
import matplotlib.pyplot as plt
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver

eta, tau, c = 0.1, 0.01, 0.5
cs, ncx, ncz, npad = 10., 25, 20, 18
hx = [(cs,ncx), (cs,npad,1.3)]
hz = [(cs,npad,-1.3), (cs,ncz), (cs,npad,1.3)]
mesh = Mesh.CylMesh([hx,1,hz], '00C')    
sigmaInf = np.ones(mesh.nC) * 0.001
actinds = mesh.gridCC[:,2]<0.
#layerind = (np.logical_and(mesh.gridCC[:,2]<0, mesh.gridCC[:,2]>-50.)) & (mesh.gridCC[:,0]<100.)
layerind = (np.logical_and(mesh.gridCC[:,2]<0, mesh.gridCC[:,2]>-50.)) & (np.logical_and(mesh.gridCC[:,0]<500, mesh.gridCC[:,0]>100.))
sigmaInf[~actinds] = 1e-8
sigmaInf[layerind] = 0.1
eta = np.zeros(mesh.nC)
eta[layerind] = 0.5
tau = np.ones(mesh.nC) * 1.
c = np.ones(mesh.nC) * 0.5
#Plot Mesh
# mesh.plotGrid()
# plt.show()


#EM calculation func
def get_em_data(sigma, eta=None, tau=None, c=None, data_type='em'):
    rxloc = np.array([[500., 0., 0.]])
    srcloc = np.array([[0., 0., 30.]])
    dt = 1.47e-3
    dt = 1.47e-3
    tpeak = 2.73e-3
    t0 = tpeak + dt
    rx_vtem = Rx.Point_dbdt(rxloc, np.logspace(np.log10(2e-5), np.log10(0.009), 51)+t0, orientation='z')
    src_vtem = EM.TDEM.Src.CircularLoop([rx_vtem], waveform=EM.TDEM.Src.VTEMWaveform(offTime=t0, peakTime=tpeak, a=3.), loc=srcloc)
    survey = Survey([src_vtem])
    if data_type == 'em':
        prb = Problem3DEM_e(mesh, sigma=sigma)
    elif data_type == 'emip':
        prb = Problem3D_e(mesh, sigmaInf=sigma, eta=eta, tau=tau, c=c)
    prb.timeSteps = [(tpeak/10, 10), ((t0-tpeak)/10, 10), (1e-06, 5), (2.5e-06, 5), (5e-06, 5), (1e-05, 10), (2e-05, 10), (4e-05, 10), (8e-05, 10), (1.6e-04, 10), (3.2e-04, 20)]
    prb.Solver = Solver
    prb.pair(survey)
    e = prb.fields(sigmaInf)
    data = survey.dpred(sigmaInf, f=e)
    # waveform
    cur = []
    for t in prb.times:
        cur.append(src_vtem.waveform.eval(t))
    cur = np.hstack(cur)
    return e, data, cur

#IP calculation func
def get_ip_data(sigma, eref, eta, tau, c):
    rxloc = np.array([[500., 0., 0.]])
    srcloc = np.array([[0., 0., 30.]])
    dt = 1.47e-3 
    rx_ip = Rx.Point_dbdt(rxloc, np.logspace(np.log10(2e-5), np.log10(0.009), 51), 'z')
    src_ip = EM.TDEM.Src.CircularLoop([rx_ip], loc=srcloc)
    dt = 1.47e-3
    tpeak = 2.73e-3
    t0 = tpeak + dt
    survey_ip = SurveyLinear([src_ip])
    t1, t2, t3 = dt, t0-0.001365, t0
    prb_ip = Problem3DIP_Linear(
        mesh, 
        sigmaInf=sigmaInf, 
        eta=eta, 
        tau=tau, 
        c=c, 
        actinds = actinds,
        tlags = [0., t1, t2, t3]
    )
    prb_ip.Solver = Solver
    prb_ip.pair(survey_ip)
    prb_ip.set_eref(eref)
    ip_approx = survey_ip.dpred([])
    return ip_approx

#Plot Model
# mesh.plotImage(np.log10(sigmaInf))
# plt.xlim(0,1000)
# plt.ylim(-200,200)
# plt.show()


e_emip, data_emip, cur = get_em_data(sigmaInf, eta, tau, c, data_type='emip')

e_em, data_em, cur = get_em_data(sigmaInf)
eref = geteref(e_em[:,0,:], mesh, option=None, tInd=20) 
ip = get_ip_data(sigmaInf, eref, eta, tau, c)
data = ip + data_em
times = np.logspace(np.log10(2e-5), np.log10(0.009), 51)

#Plot results
plt.loglog(times, -data, 'k')
plt.loglog(times, data, '--k')
plt.loglog(times, -data_emip, 'kx')
plt.loglog(times, data_emip, 'kx')
plt.loglog(times, -data_em, 'b')
plt.show()


# m = [
#     'background physical property value',
#     'layer physical property value',
#     'layer center',
#     'layer thickness'
# ]


# mesh = Mesh.TensorMesh([50, 50], x0='CC')  # 2D tensor mesh
# mapping = Maps.ParametricLayer(mesh)  # parametric layer in wholespace

# # model
# m = np.hstack(
#     np.r_[
#         1., # background value
#         2., # layer value
#         -0.1, # layer center
#         0.2 # layer thickness
#     ]
# )
# rho = mapping * m # apply the mapping


# fig, ax = plt.subplots(1, 1, figsize=(4, 6))
# mesh.plotImage(rho, ax=ax)


# plt.show()