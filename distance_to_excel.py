import numpy as np
import os
import xlsxwriter

workbook = xlsxwriter.Workbook('distance.xlsx') 
worksheet = workbook.add_worksheet() 
row = 0
column = 0

d=[]
first=-495
last=0
element=100

for i in range(element):
    distance=np.linspace(first,last,element)
    for j in range(len(distance)):
        d.append(distance[j])

    element=element-1
    last=last+5
    first=first+10

for k in range(len(d)):
    worksheet.write(row, column, d[k])
    row += 1

workbook.close()