import matplotlib.pyplot as plt
import numpy as np
import os
import xlsxwriter
from scipy import interpolate
from scipy.constants import mu_0

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
time=np.logspace(np.log10(2e-5),np.log10(2),21)
newtime=np.logspace(np.log10(2e-5),np.log10(2),200)
tt=np.logspace(np.log10(2e-5),np.log10(2),199)

rhoapp=[]


workbook = xlsxwriter.Workbook('timetest.xlsx') 
worksheet = workbook.add_worksheet() 
row = 0
column = 0

a=20
b=1
element=a
last=800
for i in range(21):
    data=np.load(os.path.join(THIS_FOLDER,'s'+str(i+1)+'.npy'))
    distance=np.linspace(40,last,element,endpoint=True)
    print(distance)
    for j in range(a):
        d=data[j+b]
        f2 = interpolate.interp1d(time, d, kind='cubic')
        ynew=f2(newtime)
        der=np.diff(ynew)
        maximum=np.argmax(abs(der))
        rhoapp.append(mu_0*distance[j]*distance[j]/tt[maximum])

    a=a-1
    b=b+1
    last=last-40
    element=element-1


for k in range(len(rhoapp)):
    worksheet.write(row, column, rhoapp[k])
    row += 1

workbook.close()