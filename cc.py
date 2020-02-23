import matplotlib.pyplot as plt
import cmath
import numpy as np

freq=np.linspace(-1, 3, 150)

eta=0.1
# tau=0.5*np.pi
tau=0.5
c=1
sigmaInf=0.01
colecole=np.empty((0,len(freq)))
phase=np.empty((0,len(freq)))

for j in range(len(freq)):
    zw=sigmaInf+sigmaInf*eta*(1/(1 + (1-eta)*(1j*2*np.pi*freq[j]*tau)**c))
    ps=zw.imag/zw.real
    phase=np.append(phase,ps)
    colecole=np.append(colecole,zw)

plt.semilogx(freq,colecole.real)
# timedom=np.real_if_close(np.fft.ifft(colecole))

# plt.plot(timedom)
# plt.xlim(0,5)
# plt.title('Cole-Cole')


# plt.plot(freq,phase)
# plt.xscale('log', basex=2)
# plt.yscale('log', basey=10)
plt.show()

# print(colecole)