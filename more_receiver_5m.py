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


rhoapp=[]

workbook = xlsxwriter.Workbook('dense_rec_eta08tau03_box.xlsx') 
worksheet = workbook.add_worksheet() 
row = 0
column = 0

a=100
b=1
element=a
last=500

def derivative(x,y):
    nx=len(x)
    dy=np.zeros(nx)
    # return (y[1:]-y[:-1])/(x[1:]-x[:-1])
    for i in range(1,nx-2):
        dy[i]=(y[i+1]-y[i-1])/(x[i+1]-x[i-1])
    
    return dy

for i in range(101):
    data=np.load(os.path.join(THIS_FOLDER,'dense_rec_eta08tau03_box_s'+str(i+1)+'.npy'))[202:303]
    distance=np.linspace(5,last,element,endpoint=True)

    for j in range(a):
        d=data[j+b]
        d=savgol_filter(d, 5, 2 , mode='nearest')
        f2 = interpolate.interp1d(time, d, kind='cubic')
        ynew=f2(newtime)

        der = 2.302*newtime*derivative(newtime,ynew)
        

        # tder=np.diff(np.log10(newtime))
        # der=np.diff(ynew)
        
        #plt.semilogx(tt,der/max(abs(der)))
        # plt.semilogx(time,d)
        # plt.semilogx(newtime,ynew,label='int')
        
        
        maximum=np.argmax(abs(der[50:]))
        
        rhoapp.append(mu_0*distance[j]*distance[j]/tt[50:][maximum])
        # print(tt[30:][maximum],rhoapp[j],distance[j])

        # fig, axs = plt.subplots(3)
        # fig.suptitle('Trimmed Data')
        # axs[0].semilogx(time,d,label='Before Interpolation')
        # axs[1].semilogx(newtime,ynew,label='After Interpolation')
        # axs[2].semilogx(tt[50:],abs(der[50:]/max(abs(der[50:]))),label=r'$\tau_{EM}$')
        # axs[0].legend(loc="upper right")
        # axs[1].legend(loc="upper right")
        # axs[2].legend(loc="upper right")

        # axs[0].set_ylabel(r'$E_x (v/m)$')
        # axs[1].set_ylabel(r'$E_x (v/m)$')
        # axs[2].set_xlabel('Time (s)')
        # axs[2].set_ylabel('Normalized pseduoimpulse response')

        # axs[0].grid()
        # axs[1].grid()
        # axs[2].grid()

        # plt.show()

    a=a-1
    b=b+1
    last=last-5
    element=element-1

# tau_em=[]
# for i in range(21):
#     data=np.load(os.path.join(THIS_FOLDER,'box_s'+str(i+1)+'.npy'))[210:231]
#     print(np.shape(data))
#     for j in range(len(data)):
#         f2 = interpolate.interp1d(time, data[j], kind='cubic')
#         ynew=f2(newtime)
#         der=np.diff(ynew)
#         maximum=np.argmax(abs(der))
#         tau_em.append(newtime[maximum])

for k in range(len(rhoapp)):
    worksheet.write(row, column, rhoapp[k])
    row += 1

workbook.close()