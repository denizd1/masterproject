import matplotlib.pyplot as plt
import numpy as np
import os
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
datafile=os.path.join(THIS_FOLDER,'blockwitheta05.npy')
datafile2=os.path.join(THIS_FOLDER,'blockwitheta01.npy')

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


phase=np.rad2deg(np.arctan2(np.asarray(imag),np.asarray(real)))
amp=np.sqrt(np.asarray(real)**2+np.asarray(imag)**2)

phase2=np.rad2deg(np.arctan2(np.asarray(imag2),np.asarray(real2)))
amp2=np.sqrt(np.asarray(real2)**2+np.asarray(imag2)**2)

plt.semilogx(freqs,phase,label='eta=0.5')
plt.semilogx(freqs,phase2,label='eta=0.1')
plt.legend()
plt.show()