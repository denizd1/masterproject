from SimPEG import Mesh, EM, Utils
import numpy as np
from simpegEMIP.TDEM.ProblemTDEMIP import Problem3D_e
from simpegEMIP.TDEM import Survey
from simpegEMIP.TDEM.Rx import Point_e
from pymatsolver import Pardiso

core_domain_x = np.r_[-300., 300.]  # extent of uniform cells in the x-direction
core_domain_z = np.r_[-200., 0.]  # extent of uniform cells in the z-direction
core_domain_y = np.r_[-250., 250.]
csx, csy, csz = 5, 15, 15

# Number of core cells in each direction
ncy = int(np.diff(core_domain_y)/csy)
ncx = int(np.diff(core_domain_x)/csx)
ncz = int(np.diff(core_domain_z)/csz)
# Number of padding cells to add in each direction
npad = 10
# Vectors of cell lengths in each direction with padding
hx = [(csx, npad, -1.5), (csx, ncx), (csx, npad, 1.5)]
hy = [(csy, npad, -1.5), (csy, ncy), (csy, npad, 1.5)]
hz = [(csz, npad, -1.5), (csz, ncz), (csz, npad, 1.5)]
# Create mesh and center it
mesh = Mesh.TensorMesh([hx, hy, hz])
#print(-mesh.hx.sum()/2)
mesh.x0 = np.r_[
    -mesh.hx.sum()/2, -mesh.hy.sum()/2., -mesh.hz[:npad+ncz].sum()
]

p0 = np.array([-100, -100, -20])
p1 = np.array([100, 100, -120])
inds = Utils.ModelBuilder.getIndicesBlock(p0, p1, mesh.gridCC)
sigma = np.ones(mesh.nC) *1e-8
sigma[mesh.gridCC[:,2]<0.] = 0.8
eta = np.ones(mesh.nC) * 0
tau = np.ones(mesh.nC) * 1e-20 
c = np.ones(mesh.nC) *1e-20

#eta[mesh.gridCC[:,2]<0.] = 0.1
#tau[mesh.gridCC[:,2]<0.] = 0.1
#c[mesh.gridCC[:,2]<0.] = 0.5
sigma[inds]=0.1
eta[inds] = 0
tau[inds] = 1e-20
c[inds] = 1e-20



x = np.linspace(-250, 250, 51)
y = np.linspace(-10, 10, 5)
rx_locations = Utils.ndgrid(x, y, np.r_[0.])
times=np.logspace(np.log10(2e-5), np.log10(10),100)

rx = Point_e(rx_locations, times, orientation='x')
pos=np.linspace(-250,250,51,endpoint=True)
for i in range(len(pos)):
	src_locations = np.array(
   		[[pos[i]-1, 0, 0],[pos[i]+1, 0, 0]]
	)
	#times = np.logspace(-3, 3, 30)
	src = EM.TDEM.Src.LineCurrent(
    		[rx], loc=src_locations,waveform=EM.TDEM.Src.StepOffWaveform(),
	)
	survey = Survey([src])
	prb_em = Problem3D_e(mesh, sigmaInf=sigma, eta=eta, tau=tau, c=c)
	# prb_em = Problem3D_e(mesh, sigma=sigma)
	prb_em.verbose = True

	#prb_em.timeSteps=[(1e-08, 10),(1e-7,10), (2.5e-06, 10), (5e-06, 10), (1e-05, 5), (2e-05, 5), (4e-05, 5), (8e-05, 5), (1.6e-04, 5), (4e-04, 5), (8e-04, 5), (1e-03, 5), (2e-03, 5), (4e-03, 5), (8e-03, 5), (1e-02, 5), (2e-02, 5), (4e-02, 5), (8e-02, 5),(5e-1,5),(8e-1,5)]


	#prb_em.timeSteps=[(1e-8,5),(1e-06, 5),(2e-05, 5),(5e-6,5),(8e-05, 5),(8e-04, 5),(2e-03, 5), (4e-03, 5),(2e-02, 5),(4e-2,5),(5e-1,5),(8e-1,5)]
	#prb_em.timeSteps = [(1e-4, 10), (1e-3, 10), (2e-2, 10), (1e1, 10), (2e2, 10), (1e3, 10), (1e5, 10)] 
	#prb_em.timeSteps = [(1e-8,10),(1e-7,10),(2.5e-6,5),(1e-4, 10), (1e-3, 10), (2e-2, 10),(5e-1,5),(8e-1,5)] old works

	prb_em.timeSteps = [(1e-8,10),(1e-7,10),(2.5e-6,10),(1e-5,10),(1e-4, 10),(4e-4,10), (1e-3, 10),(4e-3,10),(1e-2,10), (2e-2, 10),(5e-1,10),(8e-1,10)]

	prb_em.Solver = Pardiso
	prb_em.pair(survey)
	data = survey.dpred([])
	np.save('densemesh_box_s'+str(i+1)+'.npy',data.reshape((255, rx.times.size), order='F'))
