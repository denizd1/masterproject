import xlrd 
from SimPEG.Utils import plot2Ddata, ndgrid
import numpy as np

loc = r'C:\Users\dnz\Desktop\computed data\dipole-dipole\box-big_sur.xlsx'
wb = xlrd.open_workbook(loc) 
sheet = wb.sheet_by_index(0) 
data=[]
xdis=[]
ydis=[]
for i in range(sheet.nrows): 
    data.append(sheet.cell_value(i, 0))
    xdis.append(sheet.cell_value(i, 1))
    ydis.append(sheet.cell_value(i, 2))

v_max = np.max(np.abs(data))
grid=ndgrid(np.asarray(xdis),np.asarray(ydis), np.r_[0.])
cplot1 = plot2Ddata(
    grid,
    data,
    ncontour=30,
    clim=(-v_max, v_max),
    contourOpts={"cmap": "bwr"},
)