import numpy as np
from os.path import join as pjoin
import matplotlib.pyplot as plt

folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Synthetic/"
fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(8,6))

GS_file = np.load(pjoin(folder,"GS_BXZ_P.txt"))
Direct_file = np.load(pjoin(folder,"Direct_BXZ_P.txt"))
ax[0].plot(GS_file[0], c = 'b')
ax[0].plot(GS_file[1],ls=":", c = 'b')
ax[0].plot(Direct_file[0], c = 'r')
ax[0].plot(Direct_file[1],ls = ":", c = 'r')

GS_file = np.load(pjoin(folder,"GS_BXT_S.txt"))
Direct_file = np.load(pjoin(folder,"Direct_BXT_S.txt"))
ax[1].plot(GS_file[0], c = 'b')
ax[1].plot(GS_file[1],ls=":", c = 'b')
ax[1].plot(Direct_file[0], c = 'r')
ax[1].plot(Direct_file[1],ls = ":", c = 'r')

GS_file = np.load(pjoin(folder,"GS_BXR_S.txt"))
Direct_file = np.load(pjoin(folder,"Direct_BXR_S.txt"))
ax[2].plot(GS_file[0], c = 'b')
ax[2].plot(GS_file[1],ls=":", c = 'b')
ax[2].plot(Direct_file[0], c = 'r')
ax[2].plot(Direct_file[1],ls = ":", c = 'r')

GS_file = np.load(pjoin(folder,"GS_BXR_P.txt"))
Direct_file = np.load(pjoin(folder,"Direct_BXR_P.txt"))
ax[3].plot(GS_file[0], c = 'b')
ax[3].plot(GS_file[1],ls=":", c = 'b')
ax[3].plot(Direct_file[0], c = 'r')
ax[3].plot(Direct_file[1],ls = ":", c = 'r')

GS_file = np.load(pjoin(folder,"GS_BXZ_S.txt"))
Direct_file = np.load(pjoin(folder,"Direct_BXZ_S.txt"))
ax[4].plot(GS_file[0], c = 'b')
ax[4].plot(GS_file[1],ls=":", c = 'b')
ax[4].plot(Direct_file[0], c = 'r')
ax[4].plot(Direct_file[1],ls = ":", c = 'r')



plt.show()