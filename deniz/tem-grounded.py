from SimPEG import Mesh, EM, Utils
import matplotlib.pyplot as plt
import numpy as np
from simpegEMIP.TDEM.ProblemTDEMIP import Problem3D_e
from simpegEMIP.TDEM import Survey
from simpegEMIP.TDEM.Rx import Point_e
from pymatsolver import Pardiso as Solver

np.set_printoptions(threshold=np.inf)
# Cell sizes
csx, csy, csz = 50., 50., 50
# Number of core cells in each direction
ncx, ncy, ncz = 15, 15, 15
# Number of padding cells to add in each direction
npad = 4 #  need to increase this considering diffusion distance say 10
# Vectors of cell lengths in each direction with padding
hx = [(csx, npad, -1.5), (csx, ncx), (csx, npad, 1.5)]
hy = [(csy, npad, -1.5), (csy, ncy), (csy, npad, 1.5)]
hz = [(csz, npad, -1.5), (csz, ncz), (csy, npad, 1.5)]
# Create mesh and center it
mesh = Mesh.TensorMesh([hx, hy, hz], x0="CCC")


sigma = np.ones(mesh.nC) * 0.001
eta = np.ones(mesh.nC) * 0.
tau = np.ones(mesh.nC) * 1e-3
c = np.ones(mesh.nC)* 1e-8

actinds = mesh.gridCC[:,2]<0.

p0 = np.array([-100, -200, -100])
p1 = np.array([100, 200, -200])
p02 = np.array([-200, -300, -200])
p12 = np.array([200, 300, -300])

inds = Utils.ModelBuilder.getIndicesBlock(p0, p1, mesh.gridCC)
sphere=Utils.ModelBuilder.getIndicesSphere(np.array([0, 0, -500]),200,mesh.gridCC)
inds2 = Utils.ModelBuilder.getIndicesBlock(p02, p12, mesh.gridCC)

sigma[~actinds] = 1e-8

sigma[inds] = 0.1
sigma[sphere] = 0.5
sigma[inds2] = 0.9


eta[inds] = 0.1
tau[inds] = 0.7
c[inds] = 0.5

eta[sphere] = 0.6
tau[sphere] = 0.5
c[sphere] = 0.5

eta[inds2] = 0.4
tau[inds2] = 0.7
c[inds2] = 0.5



print(eta)
out = mesh.plotSlice(eta, grid=True, normal='X')
plt.colorbar(out[0])
plt.title("Chargeability (V/V)")
plt.gca().set_aspect(1)
plt.show()

x = np.linspace(-200, 200, 21)
y = np.linspace(-200, 200, 21)
rx_locations = Utils.ndgrid(x, y, np.r_[0.])

src_locations = np.array(
    [[-100, 0, 0],[100, 0, 0]]
)

# plt.plot(rx_locations[:,0], rx_locations[:,1], 'k.')
# plt.plot(src_locations[:,0], src_locations[:,1], 'r-')
# plt.gca().set_aspect(1)
# plt.show()
times = np.logspace(-3, -2, 21)

rx = Point_e(rx_locations, times, orientation='x')
src = EM.TDEM.Src.LineCurrent(
    [rx], loc=src_locations,waveform=EM.TDEM.Src.StepOffWaveform(),
)

survey = Survey([src])
prb_em = Problem3D_e(mesh, sigmaInf=sigma, eta=eta, tau=tau, c=c)
# prb_em = Problem3D_e(mesh, sigma=sigma)
prb_em.verbose = True
prb_em.timeSteps = [(1e-3, 10)]
prb_em.Solver = Solver
prb_em.pair(survey)
# m = np.r_[sigma, eta, tau, c]
# F_em = prb_em.fields(m)
data = survey.dpred([])
data1=data.reshape((21,21, rx.times.size), order='F')[15,15,:]
print(data1)


# [ 1.15890688e-05  4.28656586e-06 -3.90697719e-06 -1.31002837e-05
#  -2.34153433e-05 -3.49890305e-05 -4.79749211e-05 -4.89257687e-05
#  -4.96893610e-05 -5.05461257e-05 -5.12694833e-05 -5.17823038e-05
#  -5.23576980e-05 -5.27076669e-05 -5.30830892e-05 -5.33128569e-05
#  -5.35175583e-05 -5.36670858e-05 -5.37600635e-05 -5.38127674e-05
#  -5.38341672e-05]

# eta01dat=np.array([ 6.17189094e-05,  5.62115454e-05,  5.00321813e-05,  4.30988207e-05,
#   3.53194622e-05,  2.65908784e-05,  1.67972463e-05,  1.41694144e-05,
#   1.14071144e-05,  8.30776272e-06,  5.87358466e-06,  4.45256493e-06,
#   2.85815457e-06,  2.02280034e-06,  1.13872570e-06,  6.41871384e-07,
#   2.04412477e-07, -1.06990816e-07, -3.04132894e-07, -4.30376330e-07,
#  -5.08949692e-07])


plt.semilogx(rx.times, data1)
plt.xlabel('time')
plt.ylabel('Ex')
plt.grid(True)
plt.show()
i_time = 13


Utils.plot2Ddata(rx_locations,data.reshape((441, rx.times.size), order='F')[:,i_time])
plt.show()