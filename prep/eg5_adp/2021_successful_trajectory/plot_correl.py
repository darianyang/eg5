"""
Plot the DCC resmap from cpptraj matrix correl.
"""

import numpy as np
import matplotlib.pyplot as plt

cm = np.loadtxt("correl.dat")

plt.imshow(cm, vmin=-1, vmax=1, cmap="seismic")
#plt.imshow(cm, vmin=-1, vmax=1, cmap="bwr")
plt.colorbar()
plt.show()
