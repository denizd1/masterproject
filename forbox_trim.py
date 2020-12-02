import matplotlib.pyplot as plt
from matplotlib.widgets  import RectangleSelector
import numpy as np
import os
import xlsxwriter
from scipy import interpolate
from scipy.constants import mu_0
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import pandas as pd


THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
time=np.logspace(np.log10(2e-5),np.log10(4),200)
time2=np.logspace(np.log10(1e-5),np.log10(1),400)


rhoapp=[]

workbook = xlsxwriter.Workbook('bg_1_box_eta03_sigma10_tau500_c025.xlsx') 
worksheet = workbook.add_worksheet() 
row = 0
column = 0

a=50
b=1
element=a
last=500

def derivative(x,y):
    nx=len(x)
    dy=np.zeros(nx)
    # return (y[1:]-y[:-1])/(x[1:]-x[:-1])
    # for i in range(1,nx-1):
    #     dy[i]=(y[i+1]-y[i-1])/(x[i+1]-x[i-1])
    
    # dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
    dy[0:-1]=np.diff(y)/np.diff(x)
    dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
    return dy



#[102:153]
#[202:303]
new_path = os.path.relpath(r'.\sensitivity\bg_1_box_eta03_sigma10_tau500_c025', THIS_FOLDER)
new_path2 = os.path.relpath(r'.\sensitivity\bg_1_box_eta03_sigma10_tau5e-2_c025', THIS_FOLDER)
for i in range(51):
    data=np.load(os.path.join(new_path,'bg_1_box_eta03_sigma10_tau500_c025_s'+str(i+1)+'.npy'))
    data2=np.load(os.path.join(new_path2,'bg_1_box_eta03_sigma10_tau5e-2_c025_s'+str(i+1)+'.npy'))
    distance=np.linspace(10,last,element,endpoint=True)
    
    for j in range(a):
        
        d=data[j+b]
        d2=data2[j+b]
        d=savgol_filter(d, 65, 3, deriv=0)
        d2=savgol_filter(d2, 65, 3, deriv=0)
        der = 2.302*time2*derivative(time2,d)
        der2 = 2.302*time2*derivative(time2,d2)
        
        field=abs(der/max(abs(der)))
        field2=abs(der2/max(abs(der2)))
        
       

        # maximum=np.argmax(field)
        # maximum2=np.argmax(field2)

        # fig, axs = plt.subplots(2)
        # # fig.suptitle('Trimmed Data')
        # axs[0].semilogx(time2[:320],data[j+b][:320],label=r'$\sigma=10 S/m$')
        # axs[0].semilogx(time2[:320],data2[j+b][:320],label=r'$\sigma=12 S/m$')
        # axs[1].semilogx(time2[:320],abs(der/max(abs(der)))[:320],label=r'$\tau_{EM}$ $(\sigma=10S/m)$')
        # axs[1].semilogx(time2[:320],abs(der2/max(abs(der2)))[:320],label=r'$\tau_{EM}$ $(\sigma=12S/m)$')
        # axs[0].legend(loc="upper right")
        # axs[1].legend(loc="upper right")

        # axs[0].set_ylabel(r'$E_x (v/m)$')
        # axs[1].set_xlabel('Time (s)')
        # axs[1].set_ylabel('Normalized pseduoimpulse response')
        # axs[0].grid()
        # axs[1].grid()

        # plt.show()

        xdata=time2
        ydata=field
        fig, ax = plt.subplots()
        line, = ax.semilogx(xdata, ydata)
        point, = ax.semilogx([],[], marker="o", color="crimson")
        text = ax.text(0,0,"")
        def line_select_callback(eclick, erelease):
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

            mask= (xdata > min(x1,x2)) & (xdata < max(x1,x2)) & \
                (ydata > min(y1,y2)) & (ydata < max(y1,y2))
            xmasked = xdata[mask]
            ymasked = ydata[mask]

            if len(xmasked) > 0:
                xmax = xmasked[np.argmax(ymasked)]
                ymax = ymasked.max()
                tx = "xmax: {:.12f}, ymax {:.12f}".format(xmax,ymax)
                point.set_data([xmax],[ymax])
                text.set_text(tx)
                text.set_position((xmax,ymax))
                fig.canvas.draw_idle()
                rhoapp.append(mu_0*distance[j]*distance[j]/(6*xmax))

        rs = RectangleSelector(ax, line_select_callback,
                            drawtype='box', useblit=False, button=[1], 
                            minspanx=5, minspany=5, spancoords='pixels', 
                            interactive=True)
        plt.grid()
        plt.show()

    a=a-1
    b=b+1
    last=last-10
    element=element-1


for k in range(len(rhoapp)):
    worksheet.write(row, column, rhoapp[k])
    row += 1

workbook.close()