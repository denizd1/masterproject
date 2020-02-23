import matplotlib.pyplot as plt
import numpy as np
import os

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
datafile=os.path.join(THIS_FOLDER,'meshtest.npy')
datafile2=os.path.join(THIS_FOLDER,'meshtest3.npy')


freqs=np.logspace(-3,3,25)
real=[]
imag=[]
real2=[]
imag2=[]
data=np.load(datafile)
data2=np.load(datafile2)

for i in range(len(freqs)):
    imag.append(data[i][1])
    real.append(data[i][0])
    imag2.append(data2[i][1])
    real2.append(data2[i][0])

plt.semilogx(freqs,real2)
plt.semilogx(freqs,real)
plt.show()