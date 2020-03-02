import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.constants import mu_0
from scipy.special import erfc, erf

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))

jfile=os.path.join(THIS_FOLDER,'J05.txt')
with open(jfile, "r") as f:
    J05 = np.array([float(line) for line in f])


datafile=os.path.join(THIS_FOLDER,'freqs.txt')
with open(datafile, "r") as f:
    freqs = np.array([float(line) for line in f])


omega=freqs
j=0
mergedat=[]
real=[]
imag=[]
for i in range(10):    
    j=j+1
    jj=0
    for jj in range(4):
        jj=jj+1
        data=np.load(os.path.join(THIS_FOLDER,'t'+str(j)+'_f'+str(jj)+'.npy'))
        for k in range(len(data)):
            real.append(data[k][0])
            imag.append(data[k][1])


for rr in range(len(real)):
    mergedat.append(complex(real[rr],imag[rr]))

stepimag=[]
for e in range(len(omega)):
     r=(1/(1j*omega[e]))*np.asarray(mergedat)[e]
     stepimag.append(r)

r  = np.sqrt( 200**2. + 0**2. + 0**2.)
x=200
sigmaInf=1
time=np.logspace(-3,1,10)
ee=[]
Efield=[]
ii=1
for tt in range(len(time)):
    t=time[tt]
    theta=np.sqrt(mu_0*sigmaInf/(4*t))
    errfunc=erfc(theta*r)
    firstterm=( ( (4/np.sqrt(np.pi))*(theta**3)*(r**3) + ((6/np.sqrt(np.pi))*theta*r) )*np.exp(-theta**2*r**2)+3*errfunc )
    last=((((4/np.sqrt(np.pi))*(theta**3)*(r**3)+(2/np.sqrt(np.pi))*theta*r))*np.exp(-theta**2*r**2)+errfunc)
    ee.append((1/(4*np.pi*sigmaInf*r**3))*(firstterm*(x**2/r**2)-last))

    Efield.append(np.sqrt(2/np.pi)*np.dot(np.asarray(stepimag).imag[tt*250:250*ii],J05)/t)
    ii=ii+1

plt.semilogx(time,-1*np.asarray(Efield),label='SimPEG')
plt.semilogx(time,ee,label="W&H 2.50")
plt.xlabel('time (s)')
plt.ylabel('Ex (V/m)')
plt.legend()
plt.show()

