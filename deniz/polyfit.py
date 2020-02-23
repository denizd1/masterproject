import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.constants import mu_0

# THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
# data=os.path.join(THIS_FOLDER,'wholespace.npy')
# datafreq=np.load(data)

# freqs = np.logspace(-3,3,25)
# datdat=[]
# for j in range(len(datafreq)):
#     a=datafreq[j][0]*1/(1j*2*np.pi*freqs[j])
#     b=1j*datafreq[j][1]*1/(1j*2*np.pi*freqs[j])
#     datdat.append(a+b)

# data=np.asarray(datdat).reshape(len(datafreq),1)

# newfreqs=np.linspace(10**-3,10**3,1024)
# # step=newfreqs[1]-newfreqs[0]


# newdata=np.interp(newfreqs,freqs,datdat)
# t=np.linspace(10**-4,10**-1,50)

r  = np.sqrt( 200**2. + 0**2. + 0**2.)
x=200
freqs=np.linspace(-10**4,10**4,10000)

time=np.linspace(0.001,1,20)
Efield=[]

eta=0
c=0
tau=0
sigmaInf=1

for i in range(len(freqs)):
    sigma=sigmaInf*(1 - eta/(1 + (1j*2.*np.pi*freqs[i]*tau)**c))
    #sigma=sigmaInf*(1 - eta/(1 + (1-eta)*(1j*2.*np.pi*freqs[i]*tau)**c))
    k  = np.sqrt( -1j*2.*np.pi*freqs[i]*mu_0*sigma)
    front = 1 * 1 / (4. * np.pi * sigma * r**3) * np.exp(-1j*k*r)
    mid   = -k**2 * r**2 + 3*1j*k*r + 3
    E=front*((x**2 / r**2)*mid + (k**2 * r**2 -1j*k*r-1.))
    Efield.append(E*1/(1j*2*np.pi*freqs[i]))

result=[]

for t in range(len(time)):   
    s=0 
    for f in range(len(freqs)):
        s=s+Efield[f]*np.exp(2*np.pi*freqs[f]*1j*time[t])
    result.append(s)

result=np.asarray(result).reshape(len(time),1)

plt.plot(time,abs(result))
plt.show()