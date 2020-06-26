from SimPEG import Mesh, EM, Utils
import matplotlib.pyplot as plt
import numpy as np
from simpegEMIP.TDEM.ProblemTDEMIP import Problem3D_e
from simpegEMIP.TDEM import Survey
from simpegEMIP.TDEM.Rx import Point_e
from pymatsolver import Pardiso

#core_domain_x = np.r_[-1000., 500.]  # extent of uniform cells in the x-direction
core_domain_x = np.r_[-500.,500.]
core_domain_z = np.r_[-400., 0.]  # extent of uniform cells in the z-direction
core_domain_y = np.r_[-450.,450.]
csx, csy, csz = 10, 15, 10

# Number of core cells in each direction
ncy = int(np.diff(core_domain_y)/csy)
ncx = int(np.diff(core_domain_x)/csx)
ncz = int(np.diff(core_domain_z)/csz)
# Number of padding cells to add in each direction
npad = 10
#Vectors of cell lengths in each direction with padding
hx = [(csx, npad, -1.5), (csx, ncx), (csx, npad, 1.5)]
hy = [(csy, npad, -1.5), (csy, ncy), (csy, npad, 1.5)]
hz = [(csz, npad, -1.5), (csz, ncz), (csz, npad, 1.5)]
# Create mesh and center it
mesh = Mesh.TensorMesh([hx, hy, hz])
#print(-mesh.hx.sum()/2)
mesh.x0 = np.r_[
    -mesh.hx.sum()/2, -mesh.hy.sum()/2., -mesh.hz[:npad+ncz].sum()
]


def Circle2D(xc, r, n):
    theta = np.linspace(-np.pi, np.pi, n)
    x = r*np.cos(Utils.mkvc(theta))+xc[0]
    y = r*np.sin(Utils.mkvc(theta))+xc[1]

    return np.c_[x, y]

x0 = np.r_[-450, -450, 0.] 
x1 = np.r_[450, -450, 0.]
x2 = np.r_[450, 450, 0.]
x3 = np.r_[-450, 450, 0.]
pts0 = np.vstack((x0, x1, x2, x3, x0))
x0 = np.r_[-250, -250, -600.]
x1 = np.r_[250, -250, -600.]
x2 = np.r_[250, 250, -600.]
x3 = np.r_[-250, 250, -600.]
pts1 = np.vstack((x0, x1, x2, x3, x0))
pts = np.vstack((pts0, pts1))
inds_porphyry = Utils.ModelBuilder.PolygonInd(mesh, pts)
inds_overburden = np.logical_and(mesh.gridCC[:,2]<0., mesh.gridCC[:,2]>-30.)
np.random.seed(11)
x = np.random.randint(-500, 500, size=(100))
y = np.random.randint(-500, 500, size=(100))
z = np.random.randint(-40, 0, size=(100))

def get_random_blocks(xyz, mesh, dx=50, dy=50, dz=25):
    inds = []
    for i in range(xyz.shape[0]):
        p0 = np.r_[xyz[i,0]-dx/2., xyz[i,1]+dy/2., xyz[i,2]+dz/2.]
        p1 = np.r_[xyz[i,0]+dx/2., xyz[i,1]-dy/2., xyz[i,2]-dz/2.]
        inds.append(Utils.ModelBuilder.getIndicesBlock(p0, p1,mesh.gridCC))
    return np.hstack(inds)
    
inds_clay = get_random_blocks(np.c_[x, y, z], mesh)
inds_mineralization = Utils.ModelBuilder.getIndicesBlock(
    np.r_[-300, 100, -50], np.r_[-200, -100, -130.], mesh.gridCC
)
pts_0 = np.c_[Circle2D(np.r_[0, 0], 400, 20), np.zeros(20)]
pts_1 = np.c_[Circle2D(np.r_[0, 0], 150, 20), np.ones(20)*-600]
pts_halo = np.vstack((pts_0, pts_1))
inds_halo = Utils.ModelBuilder.PolygonInd(mesh, pts_halo)

pts_0 = np.c_[Circle2D(np.r_[0, 0], 250, 20), np.zeros(20)]
pts_1 = np.c_[Circle2D(np.r_[0, 0], 10, 20), np.ones(20)*-600]
pts_stock = np.vstack((pts_0, pts_1))
inds_stock = Utils.ModelBuilder.PolygonInd(mesh, pts_stock)

sig_background = 1./1000.
sig_overburden = 1./900.
sig_porphyry = 1./5000.
sig_stock = 1./4500.
sig_mineralization= 1./520.
sig_halo= 1./500.
sig_clay = 1./300.

# sig_air = 1e-8
# actind = mesh.gridCC[:, 2] < 0
sigma = np.ones(mesh.nC)*1e-8
sigma[mesh.gridCC[:,2]<0.] = sig_background
sigma[inds_porphyry] = sig_porphyry
sigma[inds_halo] = sig_halo
sigma[inds_stock] = sig_stock
sigma[inds_overburden] = sig_overburden
sigma[inds_clay] = sig_clay
sigma[inds_mineralization] = sig_mineralization
sig = [sig_background, sig_overburden, sig_porphyry, sig_stock, sig_halo, sig_clay, sig_mineralization]

geo = np.ones(mesh.nC) * np.nan
for i, sig_temp in enumerate(sig):
    geo[sigma==sig_temp] = i

eta = np.zeros(mesh.nC) + 1e-2
eta[geo==4] = 0.1
eta[geo==5] = 0.1
eta[geo==6] = 0.1
tau = np.ones(mesh.nC) *0.1
tau[geo==4] = 0.5
tau[geo==5] = 0.5
tau[geo==6] = 5.
c = np.ones(mesh.nC)
c[geo==4] = 0.5
c[geo==5] = 0.8
c[geo==6] = 0.5




x = np.linspace(-400, 400, 81)
y = np.linspace(-10, 10, 5)
rx_locations = Utils.ndgrid(x, y, np.r_[0.])

times=np.logspace(np.log10(2e-5), np.log10(10), 100)
rx = Point_e(rx_locations, times, orientation='x')
pos=np.linspace(-400,400,81,endpoint=True)
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

	# prb_em.timeSteps = [(1e-8,10),(1e-7,10),(2.5e-6,5),(1e-4, 10), (1e-3, 10), (2e-2, 10),(5e-1,5),(8e-1,5)]
	prb_em.timeSteps = [(1e-8,10),(1e-7,10),(2.5e-6,10),(1e-5,10),(1e-4, 10),(4e-4,10), (1e-3, 10),(4e-3,10),(1e-2,10), (2e-2, 10),(5e-1,10),(8e-1,10)]



	#prb_em.timeSteps = [(1e-8,10),(1e-7,10),(2.5e-6,5),(1e-4, 10), (1e-3, 10),(4e-3,5),(1e-2,5), (2e-2, 10),(5e-1,5),(8e-1,5)]
	#prb_em.timeSteps=[(1e-06, 5), (2.5e-06, 5), (5e-06, 5), (1e-05, 5), (2e-05, 5), (4e-05, 5), (8e-05, 5), (1.6e-04, 5), (4e-04, 5), (8e-04, 5), (1e-03, 5), (2e-03, 5), (4e-03, 5), (8e-03, 5), (1e-02, 5), (2e-02, 5), (4e-02, 5), (8e-02, 5), (1e-01, 5)]

	#prb_em.timeSteps = [(1e-4, 10), (1e-3, 10), (2e-2, 10), (1e1, 10), (2e2, 10), (1e3, 10), (1e5, 10)] 
	prb_em.Solver = Pardiso
	prb_em.pair(survey)
	data = survey.dpred([])
	np.save('por_s'+str(i+1)+'.npy',data.reshape((405, rx.times.size), order='F'))
