from SS_MTI.Read_H5 import Read_Direct_Inversion, Read_GS_h5
from SS_MTI.PostProcessing import Plot_Direct_BB
import SS_MTI.MTDecompose as _Decomp

import matplotlib.pyplot as plt
from os.path import join as pjoin
import numpy as np
from obspy.imaging.beachball import beachball

folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/invert_entire_window"

event_name = "S0173a"
depth = 60
fmin = 0.1
fmax = 0.7
misfit_name = "L2"
veloc_model = "TAYAK_BKE"
amount_of_phases = 2
az = 273.05  # S0235b: 257.57, S0183a:253

# file_name = f"Direct_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_model}.hdf5"
# Direct_File = pjoin(folder, file_name)


file_name = f"GS_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_model}.hdf5"
GS_file = pjoin(folder, file_name,)

## ============ Read Direct ========================
# (
#     depth_Direct,
#     MT_Direct,
#     DC_MT,
#     CLVD_MT,
#     misfit_L2_Direct,
#     Epsilon,
#     M0,
#     M0_DC,
#     M0_CLVD,
#     angles,
# ) = Read_Direct_Inversion(Direct_File, amount_of_phases=amount_of_phases)

# M = np.array(
#     [
#         [DC_MT[1], -DC_MT[5], DC_MT[3]],
#         [-DC_MT[5], DC_MT[2], -DC_MT[4]],
#         [DC_MT[3], -DC_MT[4], DC_MT[0]],
#     ]
# )
# from pyrocko import moment_tensor as mtm

# m = mtm.MomentTensor(
#     mnn=DC_MT[1], mee=DC_MT[2], mdd=DC_MT[0], mne=-DC_MT[5], mnd=DC_MT[3], med=-DC_MT[4]
# )

# (s1, d1, r1), (s2, d2, r2) = m.both_strike_dip_rake()
# print(s1, d1, r1)
# print(s2, d2, r2)
# print(M0_DC)


# fig = Plot_Direct_BB(
#     MT_Full=MT_Direct / M0,
#     Eps=Epsilon,
#     MT_DC=DC_MT / M0_DC,
#     M0_DC=M0_DC,
#     MT_CLVD=CLVD_MT / M0_CLVD,
#     M0_CLVD=M0_CLVD,
#     azimuths=[az, az, az],
#     inc_angles=angles,
#     phase_names=["P", "S", "pP"],
#     color="r",
#     height=19.0,
#     horizontal=True,
# )

# plt.savefig(
#     pjoin(folder, f"Test_Direct_BB_{event_name}_{depth}_{misfit_name}_{veloc_model}.svg",),
#     dpi=300,
# )
# plt.close()

## ============ Read GS ========================

depth_GS, sdr, M0_GS, misfit_L2_GS = Read_GS_h5(
    Filename=GS_file, amount_of_phases=amount_of_phases
)
Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
lowest_ind = Total_L2_GS.argsort()
Total_L2_GS.sort()
misfit_low = Total_L2_GS[:] - Total_L2_GS[0]
uncert = 0.05 * Total_L2_GS[0]
inds = np.where(misfit_low < uncert)
lowest_indices = lowest_ind[inds]

# n_lowest = int(len(Total_L2_GS) * 0.05)
# lowest_indices = Total_L2_GS.argsort()[0:n_lowest:50]
# n_lowest = 10
# lowest_indices = Total_L2_GS.argsort()[0:n_lowest]
MT = sdr[lowest_indices, :]
depth_GS = depth_GS[lowest_indices]
M0 = M0_GS[lowest_indices]

print(M0_GS[lowest_indices[0]])

plt.plot(M0_GS, np.sum(misfit_L2_GS, axis=1))
plt.plot(M0_GS[lowest_indices[0]], np.sum(misfit_L2_GS, axis=1)[lowest_indices[0]], "ro")
plt.show()
a = 1
