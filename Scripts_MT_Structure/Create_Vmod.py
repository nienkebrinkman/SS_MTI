from obspy.taup import TauPyModel as TauPyModel
from os.path import isfile, join
import numpy as np


def create_dat_file(
    src_depth: float,
    epi: float,
    baz: float,
    focal_mech: [float],
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
    :param epi: Epicentral distance (in degrees)
    :param baz: Back-azimuth (in degrees)
    :param save_path: path to save .dat file
    :param fdom: dominant frequency
    """

    bm_file = bm_file_path

    f = np.loadtxt(bm_file, skiprows=5)
    f_ud = np.flipud(f)

    radius_mars = f_ud[0][0]  # 3390 (km)

    # radius_of_planet = 3390
    km_per_deg = np.pi * (radius_mars * 1e-3) / 180.0
    dist_in_km = epi * np.pi * (radius_mars * 1e-3) / 180.0
    dist = dist_in_km

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

    moment_tensor = f"{M_tt_ins:10.4f}{M_tp_ins:10.4f}{-M_rt_ins+0:10.4f}{M_pp_ins:10.4f}{-M_rp_ins+0:10.4f}{M_rr_ins:10.4f}"
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
            print(i)
            if layer[0] == 0.0:
                continue
            depth = (radius_mars - layer[0]) * 1e-3
            if depth > 160:
                break
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
                # if vs != 0.0 and i != 0:
                dens0 = f_ud[i - 1, 1] * 1e-3
                vp0 = f_ud[i - 1, 2] * 1e-3
                vs0 = f_ud[i - 1, 3] * 1e-3
                qka0 = f_ud[i - 1, 4]
                qmu0 = f_ud[i - 1, 5]
                qs0 = qmu0
                L0 = (4 / 3) * (vs0 / vp0) ** 2
                qp0 = 1 / (L0 * (1 / qmu0) + (1 - L0) * (1 / qka0))
                if np.isnan(qp):
                    qp0 = qka0
                    qs0 = 10.0

                prev_depth = (radius_mars - f_ud[i - 1, 0]) * 1e-3
                layer_thickness = depth - prev_depth
                factor = 0.07
                layer_thickness_lim = factor * (
                    vs / fdom
                )  # layer limit should be less then 1/10 of wavelength
                if layer_thickness_lim > factor * (vs0 / fdom):
                    layer_thickness_lim = factor * (vs0 / fdom)
                import math

                n_layers = math.ceil(layer_thickness / layer_thickness_lim)
                from scipy import interpolate

                int_dens = interpolate.interp1d([prev_depth, depth], [dens0, dens])
                int_vp = interpolate.interp1d([prev_depth, depth], [vp0, vp])
                int_vs = interpolate.interp1d([prev_depth, depth], [vs0, vs])
                int_qka = interpolate.interp1d([prev_depth, depth], [qka0, qka])
                int_qmu = interpolate.interp1d([prev_depth, depth], [qmu0, qmu])
                for new_layer in np.linspace(0, layer_thickness, n_layers):
                    if new_layer == 0:
                        continue
                    depth1 = prev_depth + new_layer
                    dens1 = int_dens(depth1)
                    vp1 = int_vp(depth1)
                    vs1 = int_vs(depth1)
                    qka1 = int_qka(depth1)
                    qmu1 = int_qmu(depth1)

                    qs1 = qmu1
                    L1 = (4 / 3) * (vs1 / vp1) ** 2
                    qp1 = 1 / (L1 * (1 / qmu1) + (1 - L1) * (1 / qka1))
                    if np.isnan(qp):
                        qp1 = qka1
                        qs1 = 10.0
                    text = f"{depth1:10.4f}{vp0:10.4f}{qp0:10.4f}{vs0:10.4f}{qs0:10.4f}{dens0:10.4f}{1:10d}\n"
                    f.write(text)
                    vp0 = vp1
                    vs0 = vs1
                    qp0 = qp1
                    qs0 = qs1
                    dens0 = dens1
                    text = f"{depth1:10.4f}{vp1:10.4f}{qp1:10.4f}{vs1:10.4f}{qs1:10.4f}{dens1:10.4f}{0:10d}\n"
                    f.write(text)
            else:
                text = (
                    f"{depth:10.4f}{vp:10.4f}{qp:10.4f}{vs:10.4f}{qs:10.4f}{dens:10.4f}{1:10d}\n"
                )
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
        f.write("    0.0250     65536         0         2    0.0250  491.5200\n")

    f.close()


save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/"

bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
create_dat_file(
    src_depth=20.0,
    focal_mech=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    M0=None,
    epi=20.0,
    baz=0.0,
    save_path=save_path,
    bm_file_path=bm_file_path,
)
