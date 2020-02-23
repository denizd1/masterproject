from SimPEG import Mesh, Maps
import numpy as np
import matplotlib.pyplot as plt

def run(plotIt=True):
    mesh = Mesh.CylMesh([7, 6, 10])
    mesh.plotGrid()
    # mesh.plotGrid(slice='z')
    # mesh.plotGrid(slice='theta')

if __name__ == '__main__':
    run()
    plt.show()