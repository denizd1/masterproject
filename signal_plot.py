from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
#square
# y=[]
# x=np.linspace(0,1,500)
# for i in range(len(x)):
#     if x[i]==0:
#         y.append(0)
#     elif 0<x[i]<0.2:
#         y.append(1)
#     elif 0.2<x[i]<0.4:
#         y.append(0)
#     elif 0.4<x[i]<0.6:
#         y.append(-1)
#     elif 0.6<x[i]<0.8:
#         y.append(0)
#     elif 0.8<x[i]<1:
#         y.append(1)
#     else:
#         y.append(0)

# plt.plot(x, y)
# plt.xlabel('Time')
# plt.title('TDEM')
# plt.ylim(-2, 2)
# plt.grid()
# plt.show()

#harmonic
time = np.arange(0, 20, 0.1)
amplitude = np.sin(time)
plt.plot(time, amplitude)
plt.title('FDEM')
plt.xlabel('Time')
plt.grid()
plt.show()