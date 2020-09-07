from SS_MTI.Read_H5 import Read_Direct_Inversion
from SS_MTI.PostProcessing import Plot_Direct_BB
import SS_MTI.MTDecompose as _Decomp

import matplotlib.pyplot as plt
from os.path import join as pjoin
import numpy as np

folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Trial_8"

event_name = "S0173a"
depth = 8
fmin = 0.1
fmax = 0.7
misfit_name = "L2"
veloc_model = "TAYAK_BKE"
amount_of_phases = 5
az = 273.05  # S0235b: 257.57, S0183a:253

file_name = f"Direct_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_model}.hdf5"
Direct_File = pjoin(folder, file_name)

## ============ Read Direct ========================
(
    depth_Direct,
    MT_Direct,
    DC_MT,
    CLVD_MT,
    misfit_L2_Direct,
    Epsilon,
    M0,
    M0_DC,
    M0_CLVD,
    angles,
) = Read_Direct_Inversion(Direct_File, amount_of_phases=amount_of_phases)

M = np.array(
    [
        [DC_MT[1], -DC_MT[5], DC_MT[3]],
        [-DC_MT[5], DC_MT[2], -DC_MT[4]],
        [DC_MT[3], -DC_MT[4], DC_MT[0]],
    ]
)
from pyrocko import moment_tensor as mtm

m = mtm.MomentTensor(
    mnn=DC_MT[1], mee=DC_MT[2], mdd=DC_MT[0], mne=-DC_MT[5], mnd=DC_MT[3], med=-DC_MT[4]
)

(s1, d1, r1), (s2, d2, r2) = m.both_strike_dip_rake()
print(s1, d1, r1)
print(s2, d2, r2)
print(M0_DC)

a = 1

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
