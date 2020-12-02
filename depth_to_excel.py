import numpy as np
import os
import xlsxwriter

workbook = xlsxwriter.Workbook('depth.xlsx') 
worksheet = workbook.add_worksheet() 
row = 0
column = 0
d=[]
first=0
last=99
element=100

for i in range(element):
    depth=np.linspace(first,last,element,endpoint=True)
    for j in range(len(depth)):
        d.append(depth[j])
    element=element-1
    last=last-1

for k in range(len(d)):
    worksheet.write(row, column, d[k])
    row += 1

workbook.close()
