from SimPEG import *
from scipy.constants import mu_0
import numpy as np
import scipy.sparse as sp
from simpegEMIP.StretchedExponential import SEInvImpulseProblem, SESurvey
import matplotlib.pyplot as plt
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as PardisoSolver
from simpegEM1D import DigFilter

def ColeColePelton(f, sigmaInf, eta, tau, c):
    w = 2*np.pi*f
    return sigmaInf*(1 - eta/(1 + (1-eta)*(1j*w*tau)**c))

def ColeColeSeigel(f, sigmaInf, eta, tau, c):
    w = 2*np.pi*f
    return sigmaInf*(1 - eta/(1 + (1j*w*tau)**c))

time = np.logspace(-6,np.log10(0.01), 41)
wt, tbase, omega_int = DigFilter.setFrequency(time)
frequency = omega_int / (2*np.pi)
siginf, eta0, tau0, c0 = 1., 0.8, 0.005, 0.6
sigma = ColeColeSeigel(frequency, siginf, eta0, tau0, c0)
sigTCole = DigFilter.transFiltImpulse(sigma, wt, tbase, omega_int, time, tol=1e-12)
wires = Maps.Wires(('eta', 1), ('tau', 1), ('c', 1))
taumap = Maps.ExpMap(nP=1)*wires.tau
survey = SESurvey()
dtrue = -sigTCole
survey.dobs = dtrue
m1D = Mesh.TensorMesh([np.ones(3)])
prob = SEInvImpulseProblem(m1D, etaMap = wires.eta, tauMap = taumap, cMap=wires.c)
update_sens = Directives.UpdateSensitivityWeights()
prob.time = time
prob.pair(survey)
m0 = np.r_[eta0, np.log(tau0), c0]
perc = 0.05
dmisfitpeta = DataMisfit.l2_DataMisfit(survey)
dmisfitpeta.W = 1/(abs(survey.dobs)*perc)
reg = Regularization.Simple(m1D)
opt = Optimization.ProjectedGNCG(maxIter = 10)
invProb = InvProblem.BaseInvProblem(dmisfitpeta, reg, opt)
# Create an inversion object
target = Directives.TargetMisfit()
betaSch = Directives.BetaSchedule(coolingFactor=1, coolingRate=1)
invProb.beta = 0.
inv = Inversion.BaseInversion(invProb, directiveList=[betaSch, target])
reg.mref = 0.*m0
prob.counter = opt.counter = Utils.Counter()
opt.LSshorten = 0.5
opt.remember('xc')
opt.tolX = 1e-20
opt.tolF = 1e-20
opt.tolG = 1e-20
opt.eps = 1e-20
# opt.lower = 0.
mopt = inv.run(m0)

def ColeSEImpulse(time, eta, tau, c):    
    return eta*c/time*((time/tau)**c)*np.exp(-(time/tau)**c)
def ColeSEImpulseApprox(time, eta, tau, c):    
    m = eta*c/(tau**c)    
    peta = m*time**(c-1.) * (1. - (time/tau)**c )  
    return peta

nskip=1
fig = plt.figure()
plt.loglog(prob.time, survey.dobs, 'k', lw=1)
plt.loglog(prob.time, invProb.dpred, 'ko', alpha = 0.5)
vmin, vmax = survey.dobs.min(), survey.dobs.max()
plt.ylim(vmin, vmax)
# legend(("pred", "obs"))
print (("True      eta, tau and c = %.e, %.e and %3.1f")%(eta0, tau0, c0))
print (("Initial   eta, tau and c = %.e, %.e and %3.1f")%(eta0, tau0, c0))
print (("Estimated eta, tau and c = %.e, %.e and %3.1f")%(mopt[0], np.exp(mopt[1]), mopt[2]))

plt.xlabel("Time (s)")
plt.ylabel("$\\tilde {\eta}^I (t)$")
plt.legend(("Cole-Cole", "Stretched Exp."), fontsize = 12)
plt.grid(True)
plt.show()