import numpy as np
from SimPEG import Utils
import SimPEG.EM as EM
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

x = np.arange(-100.5, 100.5, step=1.)
y = np.r_[0]
z = x
XYZ = Utils.ndgrid(x, y, z)

Ex, Ey, Ez = EM.Analytics.FDEM.ElectricDipoleWholeSpace(
    XYZ, 
    srcLoc=np.r_[0., 0., 0.], 
    sig=1, 
    f=1.0
)
absE = np.sqrt(Ex*Ex.conj()+Ey*Ey.conj()+Ez*Ez.conj()).real

fig, ax = plt.subplots(1, 1, figsize=(6, 5))
bxplt = Ex.reshape(x.size, z.size)
bzplt = Ez.reshape(x.size, z.size)
pc = ax.pcolor(x, z, absE.reshape(x.size, z.size), norm=LogNorm())
ax.streamplot(x, z, bxplt.real, bzplt.real, color='k', density=1)
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([z.min(), z.max()])
ax.set_xlabel('x')
ax.set_ylabel('z')
cb = plt.colorbar(pc, ax=ax)
cb.set_label('|E|')
plt.show()