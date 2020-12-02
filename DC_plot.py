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


rhoapp=[]

workbook = xlsxwriter.Workbook('bg_1_box_eta03_sigma10_tau500_c025_DC.xlsx') 
worksheet = workbook.add_worksheet() 
row = 0
column = 0
rhoapp=[]
a=50
b=1
element=a
last=500
n=np.linspace(5,250,50)
new_path = os.path.relpath(r'.\sensitivity\bg_1_box_eta03_sigma10_tau500_c025', THIS_FOLDER)
for i in range(51):
    data=np.load(os.path.join(new_path,'bg_1_box_eta03_sigma10_tau500_c025_s'+str(i+1)+'.npy'))

    distance=np.linspace(10,last,element,endpoint=True)
    
    
    for j in range(a):
        d=data[j+b]
        # d=savgol_filter(d, 5, 2 , mode='nearest')
        # f2 = interpolate.interp1d(time, d, kind='cubic')
        # ynew=f2(newtime)
        rhoapp.append(abs(d[0]*2*np.pi*n[j]*(n[j]+1)*(n[j]+2)))
        
        
    a=a-1
    b=b+1
    last=last-10
    element=element-1
    np.delete(n,-1)



for k in range(len(rhoapp)):
    worksheet.write(row, column, rhoapp[k])
    row += 1

workbook.close()