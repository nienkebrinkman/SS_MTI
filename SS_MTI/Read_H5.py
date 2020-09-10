import numpy as _np
import h5py


def Read_GS_h5(Filename, amount_of_phases=5):
    """  Read Hdf5 files created by the Grid Search"""

    f = h5py.File(Filename, "r")
    misfit_L2 = f["samples"][:, -amount_of_phases:]
    nlowest = 10
    lowest_indices = misfit_L2.argsort()[0:nlowest]
    depth = f["samples"][:, 0]
    sdrs = f["samples"][:, 1:4]
    M0_corrs_total = f["samples"][:, 5]
    M0_total = f["samples"][:, 4]
    M0 = M0_total * M0_corrs_total
    return depth, sdrs, M0, misfit_L2


def Read_Direct_Inversion(Filename, amount_of_phases=5):
    f = h5py.File(Filename, "r")
    depth = f["samples"][0, 0]
    MT = f["samples"][0, 6 : 6 + 6]
    DC_MT = f["samples"][0, 6 + 6 : 6 + 2 * 6]
    CLVD_MT = f["samples"][0, 6 + 2 * 6 : 6 + 3 * 6]

    angles = f["samples"][0, 6 + 3 * 6 : -amount_of_phases]

    misfit_L2 = f["samples"][0, -amount_of_phases:]

    Cond_nr = f["samples"][0, 1]
    Epsilon = f["samples"][0, 2]
    M0 = f["samples"][0, 3]
    M0_DC = f["samples"][0, 4]
    M0_CLVD = f["samples"][0, 5]

    return depth, MT, DC_MT, CLVD_MT, misfit_L2, Epsilon, M0, M0_DC, M0_CLVD, angles, Cond_nr
