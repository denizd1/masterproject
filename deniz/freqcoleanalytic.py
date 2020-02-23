import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.constants import mu_0

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
#whole same sigma=1, block different eta (0.1,0.5)
datafile=os.path.join(THIS_FOLDER,'blockwitheta01.npy')
datafile05=os.path.join(THIS_FOLDER,'blockwitheta05.npy')
#whole same sigma=1, block different eta(0.5) and sigma(0.1)
datasig=os.path.join(THIS_FOLDER,'blocksig01eta05.npy')
#air sigma=1e-8 host and block sigma=1 block eta=0.5
datawithair=os.path.join(THIS_FOLDER,'airblocksig1eta05.npy')
#air sigma=1e-8 host sigma=1 and block sigma=1 block eta=0.5
datawairhost=os.path.join(THIS_FOLDER,'airblocksig01eta05.npy')


data=np.load(datafile)
data05=np.load(datafile05)
datasi=np.load(datasig)
datwair=np.load(datawithair)
datairhost=np.load(datawairhost)


#whole same sigma, block different eta
datreal=[]
datimag=[]
datreal05=[]
datimag05=[]
#whole same sigma, block different eta and sigma
datsigi=[]
datsigr=[]
#air sigma=1e-8 host and block sigma=1 block eta=0.5
wair_real=[]
wair_imag=[]
#air sigma=1e-8 host sigma=1 and block sigma=1 block eta=0.5
airhost_real=[]
airhost_imag=[]

for j in range(len(data)):
    datreal.append(data[j][0])
    datimag.append(data[j][1])
    datreal05.append(data05[j][0])
    datimag05.append(data05[j][1])
    datsigr.append(datasi[j][0])
    datsigi.append(datasi[j][1])
    wair_real.append(datwair[j][0])
    wair_imag.append(datwair[j][1])
    airhost_real.append(datairhost[j][0])
    airhost_imag.append(datairhost[j][1])


r  = np.sqrt( 200**2. + 0**2. + 0**2.)
x=200
freqs=np.logspace(-3,3,25)
Efield=[]
real=[]
imag=[]
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
    Efield.append(E)


for k in range(len(Efield)):
    real.append(Efield[k].real)
    imag.append(Efield[k].imag)


plt.semilogx(freqs,datreal,label='eta=0.1')
plt.semilogx(freqs,datreal05,label='eta=0.5')
plt.semilogx(freqs,datsigr,label='eta=0.5 sigma=0.1')
plt.semilogx(freqs,wair_real,label='eta=0.5 sigma=1 air=1e-8')
plt.semilogx(freqs,airhost_real,label='eta=0.5 sigma=0.1 hostsigma=1 air=1e-8')
plt.semilogx(freqs,real,label='W&H 2.40')
plt.xlabel('Frequency')
plt.ylabel('Real')
plt.legend()
plt.show()
