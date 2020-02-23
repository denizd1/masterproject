import SimPEG as simpeg
from SimPEG.EM import NSEM
import numpy as np
import matplotlib.pyplot as plt

try:
    from pymatsolver import Pardiso as Solver
except:
    from SimPEG import Solver

M = simpeg.Mesh.TensorMesh(
    [
        [(100, 9, -1.5), (100., 13), (100, 9, 1.5)],
        [(100, 9, -1.5), (100., 13), (100, 9, 1.5)],
        [(50, 10, -1.6), (50., 10), (50, 6, 2)]
    ], x0=['C', 'C', -14926.8217]
)
M.plotGrid()
plt.show()