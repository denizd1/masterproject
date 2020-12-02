import scipy.interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import pylab as pl  
# from mpl_toolkits.mplot3d import Axes3D

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))

data=np.load(os.path.join(THIS_FOLDER,'3D_box_eta03_sigma10_tau05_c025_s.npy'))
data2=np.load(os.path.join(THIS_FOLDER,'3D_box_eta03_sigma11_tau05_c025_s.npy'))
newdata=np.reshape(data,(26,51,51,200))
#print(data[102],newdata[0][2][0])
newdata2=np.reshape(data2,(26,51,51,200))
newdata=newdata-newdata2



img = newdata[:,:,:,0] # extracting the (351:467:300) array 
imgShape = np.shape(img) 
print(imgShape[2])
pl.ion()                               #added 
for i in range(imgShape[1]):           #no need for 0
    pl.cla()                           #added
    pl.imshow(img[i,:,:])
#    time.sleep(0.3) 
    pl.draw()
    pl.autoscale()
    pl.pause(0.3)                      #added   
pl.ioff()                              #added
print('Done!')



# import numpy as np
# import scipy.interpolate
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Plot slices of the data at the given coordinates
# def plot_slices(x, y, z, data, xslice, yslice, zslice, ax=None):
#     if ax is None:
#         ax = plt.figure().add_subplot(111, projection='3d')
#     # Normalize data to [0, 1] range
#     vmin, vmax = data.min(), data.max()
#     data_n = (data - vmin) / (vmax - vmin)
#     # Take slices interpolating to allow for arbitrary values
#     data_x = scipy.interpolate.interp1d(x, data, axis=0)(xslice)
#     data_y = scipy.interpolate.interp1d(y, data, axis=1)(yslice)
#     data_z = scipy.interpolate.interp1d(z, data, axis=2)(zslice)
#     print(np.shape(data_x),np.shape(data_y),np.shape(data_z))
#     # Pick color map
#     cmap = plt.cm.plasma
#     # Plot X slice
#     xs, ys, zs = data.shape
#     xplot = ax.plot_surface(xslice, y[:, np.newaxis], z[np.newaxis, :],
#                             rstride=1, cstride=1, facecolors=cmap(data_x), shade=False)
#     # Plot Y slice
#     yplot = ax.plot_surface(x[:, np.newaxis], yslice, z[np.newaxis, :],
#                             rstride=1, cstride=1, facecolors=cmap(data_y), shade=False)
#     # Plot Z slice
#     zplot = ax.plot_surface(x[:, np.newaxis], y[np.newaxis, :], np.atleast_2d(zslice),
#                             rstride=1, cstride=1, facecolors=cmap(data_z), shade=False)
#     return xplot, yplot, zplot

# np.random.seed(0)
# x = np.linspace(0, 26, 26)
# y = np.linspace(0, 51, 51)
# z = np.linspace(0, 51, 51)
# t = np.log(abs(newdata[:,:,:,199]))
# ax = plt.figure().add_subplot(111, projection='3d')
# plot_slices(x, y, z, t, 15, 26, 26, ax=ax)
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# plt.show()


# from IPython.display import HTML
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# data = newdata[:,:,:,199]
# print(np.shape(data))
# fig, ax = plt.subplots()


# im = ax.imshow(data[0,:,:])

# def init():
#     im.set_data(data[0,:,:])
#     return (im,)

# # animation function. This is called sequentially
# def animate(i):
#     data_slice = newdata[i,:,:,199]
#     im.set_data(data_slice)
#     return (im,)

# # call the animator. blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=24, interval=20)
# plt.show()
# #HTML(anim.to_html5_video())