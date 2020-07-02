import matplotlib.pyplot as plt
import numpy as np
import os
import xlsxwriter
from scipy import interpolate
from scipy.constants import mu_0
from scipy.signal import savgol_filter

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
time=np.logspace(np.log10(2e-5),np.log10(10),100)
newtime=np.logspace(np.log10(2e-5),np.log10(10),200)
tt=np.logspace(np.log10(2e-5),np.log10(10),199)


a=50
b=1
element=a
last=500
for i in range(51):
    data1=np.load(os.path.join(THIS_FOLDER,'box_s'+str(i+1)+'.npy'))[102:153]
    data2=np.load(os.path.join(THIS_FOLDER,'eta04tau03_box_s'+str(i+1)+'.npy'))[102:153]
    data3=np.load(os.path.join(THIS_FOLDER,'eta04_box_s'+str(i+1)+'.npy'))[102:153]
    data4=np.load(os.path.join(THIS_FOLDER,'eta04tau07_box_s'+str(i+1)+'.npy'))[102:153]
    data5=np.load(os.path.join(THIS_FOLDER,'eta04tau10_box_s'+str(i+1)+'.npy'))[102:153]
    data6=np.load(os.path.join(THIS_FOLDER,'eta04tau100_box_s'+str(i+1)+'.npy'))[102:153]
    distance=np.linspace(10,last,element,endpoint=True)


    for j in range(a):
        d1=data1[j+b]
        d2=data2[j+b]
        d3=data3[j+b]
        d4=data4[j+b]
        d5=data5[j+b]
        d6=data6[j+b]
        d1=savgol_filter(d1, 5, 2 , mode='nearest')
        d2=savgol_filter(d2, 5, 2 , mode='nearest')
        d3=savgol_filter(d3, 5, 2 , mode='nearest')
        d4=savgol_filter(d4, 5, 2 , mode='nearest')
        d5=savgol_filter(d5, 5, 2 , mode='nearest')
        d6=savgol_filter(d6, 5, 2 , mode='nearest')
        f1 = interpolate.interp1d(time, d1, kind='quadratic')
        f2 = interpolate.interp1d(time, d2, kind='quadratic')
        f3 = interpolate.interp1d(time, d3, kind='quadratic')
        f4 = interpolate.interp1d(time, d4, kind='quadratic')
        f5 = interpolate.interp1d(time, d5, kind='quadratic')
        f6 = interpolate.interp1d(time, d6, kind='quadratic')

        ynew1=f1(newtime)
        ynew2=f2(newtime)
        ynew3=f3(newtime)
        ynew4=f4(newtime)
        ynew5=f5(newtime)
        ynew6=f6(newtime)

        tder=np.diff(np.log10(newtime))
        
        der1=np.diff(ynew1)/tder
        der2=np.diff(ynew2)/tder
        der3=np.diff(ynew3)/tder
        der4=np.diff(ynew4)/tder
        der5=np.diff(ynew5)/tder
        der6=np.diff(ynew6)/tder
        

        fig, axs = plt.subplots(3)
        fig.suptitle('Trimmed Data')
        #fig.suptitle('Original Data')
        axs[0].semilogx(time,d1,label='Before Interpolation (No IP)')
        axs[0].semilogx(time,d2,label=r'Before Interpolation ($\eta$=0.4 $\tau$=0.3)',linestyle='--')
        axs[0].semilogx(time,d3,label=r'Before Interpolation ($\eta$=0.4 $\tau$=0.5)',linestyle='-.')
        axs[0].semilogx(time,d4,label=r'Before Interpolation ($\eta$=0.4 $\tau$=0.7)',linestyle='-.')
        axs[0].semilogx(time,d5,label=r'Before Interpolation ($\eta$=0.4 $\tau$=10)',linestyle='-.')
        axs[0].semilogx(time,d6,label=r'Before Interpolation ($\eta$=0.4 $\tau$=100)',linestyle='-.')

        axs[1].semilogx(newtime,ynew1,label='After Interpolation (No IP)')
        axs[1].semilogx(newtime,ynew2,label=r'After Interpolation ($\eta$=0.4 $\tau$=0.3)',linestyle='--')
        axs[1].semilogx(newtime,ynew3,label=r'After Interpolation ($\eta$=0.4 $\tau$=0.5)',linestyle='-.')
        axs[1].semilogx(newtime,ynew4,label=r'After Interpolation ($\eta$=0.4 $\tau$=0.7)',linestyle='-.')
        axs[1].semilogx(newtime,ynew5,label=r'After Interpolation ($\eta$=0.4 $\tau$=10)',linestyle='-.')
        axs[1].semilogx(newtime,ynew6,label=r'After Interpolation ($\eta$=0.4 $\tau$=100)',linestyle='-.')


        axs[2].semilogx(tt[30:],abs(der1[30:]/max(abs(der1[30:]))),label=r'$\tau_{EM}$ (No IP)')
        axs[2].semilogx(tt[30:],abs(der2[30:]/max(abs(der2[30:]))),label=r'$\tau_{EM}$ ($\eta$=0.4 $\tau$=0.3)',linestyle='--')
        axs[2].semilogx(tt[30:],abs(der3[30:]/max(abs(der3[30:]))),label=r'$\tau_{EM}$ ($\eta$=0.4 $\tau$=0.5)',linestyle='-.')
        axs[2].semilogx(tt[30:],abs(der4[30:]/max(abs(der4[30:]))),label=r'$\tau_{EM}$ ($\eta$=0.4 $\tau$=0.7)',linestyle='-.')
        axs[2].semilogx(tt[30:],abs(der5[30:]/max(abs(der5[30:]))),label=r'$\tau_{EM}$ ($\eta$=0.4 $\tau$=10)',linestyle='-.')
        axs[2].semilogx(tt[30:],abs(der6[30:]/max(abs(der6[30:]))),label=r'$\tau_{EM}$ ($\eta$=0.4 $\tau$=100)',linestyle='-.')
        
        axs[0].legend(loc="upper right")
        axs[1].legend(loc="upper right")
        axs[2].legend(loc="upper right")

        axs[0].set_ylabel(r'$E_x (v/m)$')
        axs[1].set_ylabel(r'$E_x (v/m)$')
        axs[2].set_xlabel('Time (s)')
        axs[2].set_ylabel('Normalized pseduoimpulse response')

        axs[0].grid()
        axs[1].grid()
        axs[2].grid()

        plt.show()

    a=a-1
    b=b+1
    last=last-10
    element=element-1

