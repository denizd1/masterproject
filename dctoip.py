"""
2.5D DC-IP inversion of Dipole Dipole array with Topography
===========================================================

This is an example for 2.5D DC-IP Inversion.
For DC inversion, a resisistivity model (Ohm-m) is generated having conductive
and resistive cylinders; they are respectively located right and left sides
of the subsurface.
For IP inversion, a chargeability model (V/V) is generated having a chargeable
cyinder located at the center.
Default `survey_type` is dipole-dipole, but this can be changed.
to 'pole-dipole', 'dipole-pole', and 'pole-pole'.
By running DC and IP simulations synthetic DC and IP data are generated,
respectively. Following two-stage approach (Oldenburg et al, 1999),
first DC data is inverted to recover a resistivity model. Then by using
the obtained resistivity model, sensitivity function is formed and used for
subsequent IP inversion to recover a chargeability model.
"""

from SimPEG import DC, IP
from SimPEG import Maps, Utils
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from pylab import hist
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver

survey_type='pole-dipole'
np.random.seed(1)
# Initiate I/O class for DC
IO = DC.IO()
# Obtain ABMN locations

xmin, xmax = 0., 200.
ymin, ymax = 0., 0.
zmin, zmax = 0, 0
endl = np.array([[xmin, ymin, zmin], [xmax, ymax, zmax]])
# Generate DC survey object
survey_dc = DC.Utils.gen_DCIPsurvey(endl, survey_type=survey_type, dim=2,
                                    a=10, b=10, n=10)
survey_dc.getABMN_locations()
survey_dc = IO.from_ambn_locations_to_survey(
    survey_dc.a_locations, survey_dc.b_locations,
    survey_dc.m_locations, survey_dc.n_locations,
    survey_type, data_dc_type='volt', data_ip_type='volt'
)

# Obtain 2D TensorMesh
mesh, actind = IO.set_mesh()
topo, mesh1D = DC.Utils.genTopography(mesh, -10, 0, its=100)
actind = Utils.surface2ind_topo(mesh, np.c_[mesh1D.vectorCCx, topo])
survey_dc.drapeTopo(mesh, actind, option="top")

# Build conductivity and chargeability model
blk_inds_c = Utils.ModelBuilder.getIndicesSphere(
    np.r_[60., -25.], 12.5, mesh.gridCC
)
blk_inds_r = Utils.ModelBuilder.getIndicesSphere(
    np.r_[140., -25.], 12.5, mesh.gridCC
)
blk_inds_charg = Utils.ModelBuilder.getIndicesSphere(
    np.r_[100., -25], 12.5, mesh.gridCC
)
sigma = np.ones(mesh.nC)*1./100.
sigma[blk_inds_c] = 1./10.
sigma[blk_inds_r] = 1./1000.
sigma[~actind] = 1./1e8
rho = 1./sigma
charg = np.zeros(mesh.nC)
charg[blk_inds_charg] = 0.1


# Use Exponential Map: m = log(rho)
actmap = Maps.InjectActiveCells(
    mesh, indActive=actind, valInactive=np.log(1e8)
)
mapping = Maps.ExpMap(mesh) * actmap

# Generate mtrue_dc for resistivity
mtrue_dc = np.log(rho[actind])

# Generate 2.5D DC problem
# "N" means potential is defined at nodes
prb = DC.Problem2D_N(
    mesh, rhoMap=mapping, storeJ=True,
    Solver=Solver
)
# Pair problem with survey
try:
    prb.pair(survey_dc)
except:
    survey_dc.unpair()
    prb.pair(survey_dc)

# Make synthetic DC data with 5% Gaussian noise
dtrue_dc = survey_dc.makeSyntheticData(mtrue_dc, std=0.05, force=True)
IO.data_dc = dtrue_dc

# Generate mtrue_ip for chargability
mtrue_ip = charg[actind]
# Generate 2.5D DC problem
# "N" means potential is defined at nodes
prb_ip = IP.Problem2D_N(
    mesh, etaMap=actmap, storeJ=True, rho=rho,
    Solver=Solver
)
survey_ip = IP.from_dc_to_ip_survey(survey_dc, dim="2.5D")
prb_ip.pair(survey_ip)

dtrue_ip = survey_ip.makeSyntheticData(mtrue_ip, std=0.05)

IO.data_ip = dtrue_ip





# Set initial model based upon histogram
m0_dc = np.ones(actmap.nP)*np.log(100.)
# Set uncertainty
# floor
eps_dc = 10**(-3.2)
# percentage
std_dc = 0.05

mopt_dc, pred_dc = DC.run_inversion(
    m0_dc, survey_dc, actind, mesh, std_dc, eps_dc,
    beta0_ratio=1e0,
    use_sensitivity_weight=True
    )

# Convert obtained inversion model to resistivity
# rho = M(m), where M(.) is a mapping

rho_est = mapping*mopt_dc
rho_est[~actind] = np.nan
rho_true = rho.copy()
rho_true[~actind] = np.nan



# Set initial model based upon histogram
m0_ip = np.ones(actmap.nP)*1e-10
# Set uncertainty
# floor
eps_ip = 10**(-4)
# percentage
std_ip = 0.05
# Clean sensitivity function formed with true resistivity
prb_ip._Jmatrix = None
# Input obtained resistivity to form sensitivity
prb_ip.rho = mapping*mopt_dc
mopt_ip, _ = IP.run_inversion(
    m0_ip, survey_ip, actind, mesh, std_ip, eps_ip,
    upper=np.Inf, lower=0.,
    beta0_ratio=1e0,
    use_sensitivity_weight=True
)

# Convert obtained inversion model to chargeability
# charg = M(m), where M(.) is a mapping for cells below topography

charg_est = actmap*mopt_ip
charg_est[~actind] = np.nan
charg_true = charg.copy()
charg_true[~actind] = np.nan