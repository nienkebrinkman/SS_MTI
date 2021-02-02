from obspy.taup import TauPyModel as TauPyModel
from obspy.geodetics import kilometers2degrees as kd
from os.path import isfile, join
import numpy as np
import re


def read_depth_from_dat(dat_folder: str):

    """
    Reads out the depth of the event specified in .dat file
    :param dat_folder: folder where .dat file is located
    :returns depth: the depth of the event
    """

    with open(join(dat_folder, "crfl.dat"), "r+") as f:
        data = f.readlines()
        f.close()
    return np.array(re.findall("\d+\.\d+", data[-9]), dtype=float)[2]


def read_epi_from_dat(dat_folder: str, radius: float = 3389.5):
    """
    Reads out the depth of the event specified in .dat file
    :param dat_folder: folder where .dat file is located
    :returns depth: the depth of the event
    """

    with open(join(dat_folder, "crfl.dat"), "r+") as f:
        data = f.readlines()
        f.close()
    return kd(np.array(re.findall("\d+\.\d+", data[-6]), dtype=float)[0], radius)


def update_dat_file(
    dat_folder: str,
    m: np.array,
    vpvs: bool,
    depth: bool,
    produce_tvel: bool = True,
    tvel_name: str = "Test",
):
    """ 
    This function will update an existing dat file with only changes 
    in the moment tensor and the Vp and Vs
    :param dat_folder: folder to your .dat file
    :param m: vector including model parameters (fmech(6), structure(x))
    :param vpvs: if true, vpvs updates will be done (starts from depth layer zero)
    :param depth: if true, depth layer updates will be done 
                  (if only 1 depth given: MOHO will change)
    :param produce_tvel: if true, a .tvel file will be produced
                         based on the updated .dat file
    :param tvel_name: name of the .tvel file
    """
    fmech = m[:6]

    with open(join(dat_folder, "crfl.dat"), "r+") as f:
        data = f.readlines()
        skiprows = 3  # Always need to skip first 3 lines
        """ Updating the moment tensor in .dat file"""
        fmech_update = f"{fmech[0]:10.4f}{fmech[1]:10.4f}{fmech[2]:10.4f}{fmech[3]:10.4f}{fmech[4]:10.4f}{fmech[5]:10.4f}\n"
        data[-8] = fmech_update
        """ Updating the structural parameters in .dat file """
        if vpvs and depth == False:
            print("vpvs are changed in dat file starting from depth 0")
            n_params = int((len(m) - 6) / 2)
            vp = m[6 : 6 + n_params]
            vs = m[6 + n_params : 6 + 2 * n_params]
            for i in range(n_params):
                line = data[skiprows + i * 2]
                """ search for and create floats from the string line """
                flt = np.array(re.findall("\d+\.\d+", line), dtype=float)
                """ replace vp and vs """
                text = f"{flt[0]:10.4f}{vp[i]:10.4f}{flt[2]:10.4f}{vs[i]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                data[skiprows + i * 2] = text
                if i != n_params - 1:
                    line = data[skiprows + i * 2 + 1]
                    """ search for and create floats from the string line """
                    flt1 = np.array(re.findall("\d+\.\d+", line), dtype=float)
                    """ replace depth """
                    text = f"{flt1[0]:10.4f}{vp[i]:10.4f}{flt[2]:10.4f}{vs[i]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                    data[skiprows + i * 2 + 1] = text
        elif depth and vpvs == False:
            n_params = int((len(m) - 6))
            depth = m[6 : 6 + n_params]
            if n_params == 1:
                print("depth of MOHO (from TAYAK) will be changed")
                flt = np.array(re.findall("\d+\.\d+", data[9]), dtype=float)
                data[
                    9
                ] = f"{depth[0]:10.4f}{flt[1]:10.4f}{flt[2]:10.4f}{flt[3]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                flt = np.array(re.findall("\d+\.\d+", data[8]), dtype=float)
                data[
                    8
                ] = f"{depth[0]:10.4f}{flt[1]:10.4f}{flt[2]:10.4f}{flt[3]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
            else:
                print("depths are changed in dat file starting from depth 0")
                for i in range(n_params):
                    line = data[skiprows + i * 2]
                    """ search for and create floats from the string line """
                    flt = np.array(re.findall("\d+\.\d+", line), dtype=float)
                    """ replace vp and vs """
                    text = f"{depth[i]:10.4f}{flt[1]:10.4f}{flt[2]:10.4f}{flt[3]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                    data[skiprows + i * 2] = text
                    if i != n_params - 1:
                        line = data[skiprows + i * 2 + 1]
                        """ search for and create floats from the string line """
                        flt1 = np.array(re.findall("\d+\.\d+", line), dtype=float)
                        """ replace depth """
                        text = f"{depth[i+1]:10.4f}{flt[1]:10.4f}{flt[2]:10.4f}{flt[3]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                        data[skiprows + i * 2 + 1] = text

        elif depth and vpvs:
            n_params = int((len(m) - 6) / 3)
            depth = m[6 : 6 + n_params]
            vp = m[6 + n_params : 6 + 2 * n_params]
            vs = m[6 + 2 * n_params : 6 + 3 * n_params]
            for i in range(n_params):
                line = data[skiprows + i * 2]
                """ search for and create floats from the string line """
                flt = np.array(re.findall("\d+\.\d+", line), dtype=float)
                """ replace vp and vs """
                text = f"{depth[i]:10.4f}{vp[i]:10.4f}{flt[2]:10.4f}{vs[i]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                data[skiprows + i * 2] = text
                if i != n_params - 1:
                    line = data[skiprows + i * 2 + 1]
                    """ search for and create floats from the string line """
                    flt1 = np.array(re.findall("\d+\.\d+", line), dtype=float)
                    """ replace depth """
                    text = f"{depth[i+1]:10.4f}{vp[i]:10.4f}{flt[2]:10.4f}{vs[i]:10.4f}{flt[4]:10.4f}{flt[5]:10.4f}{1:10d}\n"
                    data[skiprows + i * 2 + 1] = text
        else:
            raise ValueError("vpvs either True or False & depth either True or False")

        f.close()
    with open(join(dat_folder, "crfl.dat"), "w") as f:
        f.write("".join(data))
        f.close()
    """ automatically create .tvel file """
    depth = np.zeros(len(data[3:-18]))
    vp = np.zeros(len(data[3:-18]))
    vs = np.zeros(len(data[3:-18]))
    dens = np.zeros(len(data[3:-18]))
    for i, line in enumerate(data[3:-18]):
        """ search for and create floats from the string line """
        flt = np.array(re.findall("\d+\.\d+", line), dtype=float)
        depth[i] = flt[0]
        vp[i] = flt[1]
        vs[i] = flt[3]
        dens[i] = flt[5]
    create_tvel_file(depth, vp, vs, dens, dat_folder, tvel_name)


def create_tvel_file(
    depth: np.array,
    vp: np.array,
    vs: np.array,
    dens: np.array,
    save_folder: str,
    name: str = "Test",
):
    """
    Creating a .tvel file that can be used to compute a .npz file
    The vectors depth, vp, vs and dens MUST be same length
    :param depth: vector with layer depths
    :param vp: vector with P-velocity values at each depth layer
    :param vs: vector with S-velocity values at each depth layer
    :param dens: vector with density values at each depth layer
    :param save_folder: folder where .tvel file will be saved
    :param name: name of .tvel file
    """

    assert (
        len(depth) == len(vp) and len(depth) == len(vs) and len(depth) == len(dens)
    ), "All arrays (depth, vp, vs and dens) should be of same length"

    """ combining all the data vector """
    data = np.vstack((np.vstack((np.vstack((depth, vp)), vs)), dens)).T

    with open(join(save_folder, f"{name}.tvel"), "w") as f:
        f.write("# Input file for TauP\n")
        f.write("NAME         TAYAK_BKE\n")
        for line in data:
            f.write(f"{line[0]:8.2f}{line[1]:8.3f}{line[2]:8.3f}{line[3]:8.3f}\n")
        f.write(
            """ 1596.98   4.986   0.000   5.855
 1853.05   5.150   0.000   6.025
 2109.13   5.284   0.000   6.166
 2365.20   5.393   0.000   6.280
 2621.27   5.475   0.000   6.368
 2877.35   5.534   0.000   6.430
 3133.42   5.569   0.000   6.467
 3389.50   5.569   0.000   6.467"""
        )
        f.close()


def create_dat_file(
    src_depth: float,
    epi_in_km: float,
    baz: float,
    focal_mech: [float],
    dt: float,
    M0: float = None,
    save_path: str = "./",
    bm_file_path: str = "./",
    fdom: str = 1.000,
):
    """ 
    This function creates a .dat file that is used for the reflectivity code of Fuchs&Muller
    :paran src_depth: source depth
    :param focal_mech: strike,dip,rake or m_rr, m_tt, m_pp, m_rt, m_rp, m_tp
    :param M0: scalar moment, only necessesary when focal_mech strike,dip,rake
    :param epi_in_km : epi_in_km central distance (in degrees)
    :param baz: Back-azimuth (in degrees)
    :param save_path: path to save .dat file
    :param fdom: dominant frequency
    """

    bm_file = bm_file_path

    f = np.loadtxt(bm_file, skiprows=5)
    f_ud = np.flipud(f)

    radius_mars = 3389.5 * 1e3  # f_ud[0][0]  # 3390 (km)

    # radius_of_planet = 3390
    # km_per_deg = np.pi * (radius_mars * 1e-3) / 180.0
    # dist_in_km = epi_in_km  * np.pi * (radius_mars * 1e-3) / 180.0
    dist = epi_in_km

    if baz < 0:
        baz *= -1
    rec_az = baz
    rec_z = 0.0

    src_x = 0.0
    src_y = 0.0
    src_z = src_depth
    or_time = 0.0
    s_strength = 1.0

    assert (M0 is None and len(focal_mech) == 6) or (M0 is not None and len(focal_mech) == 3), (
        "focal_mech length is incorrect. "
        "If you specify M0, focal_mech is [strike,dip,rake]. "
        "Otherwise focal_mech is [m_rr, m_tt, m_pp, m_rt, m_rp, m_tp]"
    )

    for i in range(len(focal_mech)):
        focal_mech[i] += 0

    M_tt_ins = focal_mech[1]
    M_pp_ins = focal_mech[2]
    M_rr_ins = focal_mech[0]
    M_rp_ins = focal_mech[4]
    M_rt_ins = focal_mech[3]
    M_tp_ins = focal_mech[5]

    moment_tensor = f"{M_tt_ins:10.4f}{-M_tp_ins+0:10.4f}{M_rt_ins:10.4f}{M_pp_ins:10.4f}{-M_rp_ins+0:10.4f}{M_rr_ins:10.4f}"
    # moment_tensor = f"{M_tt_ins:10.4f}{M_tp_ins:10.4f}{M_rt_ins:10.4f}{M_pp_ins:10.4f}{M_rp_ins:10.4f}{M_rr_ins:10.4f}"

    # model = TauPyModel(taup_path)
    # model_layers = model.model.s_mod.v_mod.layers

    with open(join(save_path, "crfl.dat"), "w") as f:
        f.write("Test name\n")
        f.write(" 0 0 0 0 0   0 0 1 1 1   2 1 0 0 1   0 1 2 0 1   1\n")
        f.write("    5    1    0    1    1\n")

        # Get the indices of the velocity model with blocky description
        indices = np.setdiff1d(
            np.arange(len(f_ud[:, 0])), np.unique(f_ud[:, 0], return_index=True)[1]
        )
        indices1 = indices - 1
        inds = np.sort(np.hstack((0, np.hstack((indices1, indices)))))

        for i, layer in enumerate(f_ud):
            if layer[0] == 0.0:
                continue
            depth = (radius_mars - layer[0]) * 1e-3
            dens = layer[1] * 1e-3
            vp = layer[2] * 1e-3
            vs = layer[3] * 1e-3
            qka = layer[4]  # qka
            qmu = layer[5]  # qmu
            vph = layer[6]
            vsh = layer[7]
            eta = layer[8]

            qs = qmu
            L = (4 / 3) * (vs / vp) ** 2
            qp = 1 / (L * (1 / qmu) + (1 - L) * (1 / qka))
            if np.isnan(qp):
                qp = qka
                qs = 10.0

            # Check if part of velocity model is part of the gradient:
            if i not in inds and vs != 0.0:
                # prev_depth = (radius_mars - f_ud[i - 1, 0]) * 1e-3
                # layer_thickness = depth - prev_depth
                # factor = 0.07
                # layer_thickness_lim = factor * (
                #     vs / fdom
                # )  # layer limit should be less then 1/10 of wavelength
                # vs0 = f_ud[i - 1, 3] * 1e-3
                # if layer_thickness_lim > factor * (vs0 / fdom):
                #     layer_thickness_lim = factor * (vs0 / fdom)
                # import math

                # n_layers = math.ceil(layer_thickness / layer_thickness_lim)
                n_layers = 1
            else:
                n_layers = 1
            text = f"{depth:10.4f}{vp:10.4f}{qp:10.4f}{vs:10.4f}{qs:10.4f}{dens:10.4f}{n_layers:10d}\n"
            f.write(text)
        f.write("\n")
        f.write(f"{rec_z:10.4f}\n")
        f.write(f"{src_x:10.4f}{src_y:10.4f}{src_z:10.4f}{or_time:10.4f}{s_strength:10.4f}\n")
        f.write(f"{moment_tensor}\n")
        f.write(f"{dist:10.4f}{dist:10.4f}{0.:10.4f}{rec_az:10.4f}{1:10d}\n")
        f.write(f"{dist:10.4f}\n")
        f.write(f"{rec_az:10.4f}\n")
        f.write(f"{12.:10.4f}   {-300.:10.4f}\n")
        f.write("    3.0000    3.5000   23.5000   25.0000      650\n")
        f.write(f"    0.0100    0.0133{fdom:10.4f}    1.0300    0.0000\n")
        # f.write("    0.2420     32768         0         2    0.2420  245.7600\n")
        npts = 32768
        t_sigma = 0.3 * dt * npts
        f.write(f"{dt:10.4f}{npts:10d}{0:10d}{2:10d}{dt:10.4f}{t_sigma:10.4f}\n")

    f.close()


def plot_dat_file(dat_paths: [str]):
    """ 
    This function plots the list of .dat files for Vp, Vs and Density over depth
    :paran dat_path: List of file_paths to your dat files
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, sharey="all", sharex="col", figsize=(8, 6))
    for i, dat_path in enumerate(dat_paths):
        if i == i:
            skipfoot = 11 + 9
        else:
            skipfoot = 11
        dat_file = pd.read_csv(
            dat_path,
            skiprows=3,
            skipfooter=skipfoot,
            header=None,
            delim_whitespace=True,
            engine="python",
        )
        depth = dat_file.values[:, 0]
        vp = dat_file.values[:, 1]
        vs = dat_file.values[:, 3]
        dens = dat_file.values[:, 5]

        ax[0].plot(vp, depth, label=f"nr {i}")

        ax[1].plot(vs, depth)
        ax[2].plot(dens, depth)
    ax[0].set_ylim(ax[0].get_ylim()[::-1])
    ax[0].legend()
    plt.show()


# save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_5/"

# bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
# create_dat_file(
#     src_depth=20.0,
#     focal_mech=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
#     M0=None,
#     epi_in_km =20.0,
#     baz=0.0,
#     save_path=save_path,
#     bm_file_path=bm_file_path,
# )

# # plot_files = [
# #     "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_2/crfl.dat",
# #     "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_3/crfl.dat",
# #     "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_4/crfl.dat",
# # ]
# plot_files = [
#     "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_5/crfl.dat",
# ]
# plot_dat_file(plot_files)
