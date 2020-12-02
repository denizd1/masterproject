import numpy as np
from scipy.special import erf
from scipy.constants import mu_0
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors


def derivative(x,y):
    nx=len(x)
    dy=np.zeros(nx)
    dy[0:-1]=np.diff(y)/np.diff(x)
    dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
    return dy

sigma=10
time=np.logspace(np.log10(1e-5),np.log10(1),500)
distance=np.linspace(10,50,51)
y=0
Ex=np.empty([len(distance), len(time)])
der_Ex=np.empty([len(distance), len(time)])
for i in range(len(distance)):
    for j in range(len(time)):
        r=np.sqrt(distance[i]**2+y**2)
        theta=np.sqrt(mu_0*sigma/(4*time[j]))
        first_term=1/(2*np.pi*sigma*r**3)
        second_term=(erf(theta*r)-(2/np.sqrt(np.pi))*theta*r*np.exp((-theta**2)*(r**2)))
        Ex[i,j]=-first_term*second_term

for i in range(len(distance)):
    der_Ex[i,:]=derivative(time,Ex[i,:])

der_Ex=2.302*time*der_Ex

for j in range(len(distance)):
    der_Ex[j,:]=abs(der_Ex[j,:]/max(abs(der_Ex[j,:])))


scl=0.4
dscl=1.5
x=np.reshape(np.linspace(1,15,51),[51,1])
plt.xlim(1e-4,1e-1)
plt.semilogx(time,(der_Ex[10:,:]+scl*x[10:]).T)

ax=plt.gca()
# label_format = '{:,.0f}'
# ax.set_yticks(ax.get_yticks().tolist())
# ax.set_yticklabels([label_format.format(x) for x in ax.get_yticks().tolist()])

ax.set_yticklabels([18,14,23,32,41,51])
plt.xlabel('Time (s)')
plt.ylabel('# of reciever')
plt.grid()
plt.show()