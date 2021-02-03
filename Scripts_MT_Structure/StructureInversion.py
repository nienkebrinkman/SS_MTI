"""
Invert focal mechanism and structure simultaneously
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from os.path import isfile
from os import listdir as lsdir
from os import makedirs
from obspy.taup import TauPyModel
import obspy
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
import toml
from obspy import UTCDateTime as utct
from scipy.optimize import approx_fprime as af

import SS_MTI
import Create_Vmod
import subprocess
from SS_MTI import PreProcess as _PreProcess
from SS_MTI import GreensFunctions as _GreensFunctions
from SS_MTI import Gradient as _Gradient


def define_arguments():
    helptext = "Simultaneous source-structure inversion"
    parser = argparse.ArgumentParser(description=helptext)

    helptext = "Input toml file"
    parser.add_argument("input_file", help=helptext)
    return parser.parse_args()


if __name__ == "__main__":
    args = define_arguments()
    print(f"Inversion based on input file: {args.input_file}")
    input_file = toml.load(args.input_file, _dict=dict)

    """ 1.Read in the observed data: """

    """ 1.1 Read the inventory and catalog file (the once that contain info about the marsquakes) """
    inv = SS_MTI.DataGetter.read_inv(
        inv_path=input_file["DATA"]["inventory_filepath"]
    )  # Inventory file
    cat = SS_MTI.DataGetter.read_cat(
        cat_path=input_file["DATA"]["catalog_filepath"]
    )  # Catalog file

    """ 1.2 Get the data into a list of obspy.Event objects """
    save_folder = input_file["SAVE"]["save_folder"]
    events = SS_MTI.DataGetter.read_events_from_cat(
        event_params=input_file["EVENTS"],
        cat=cat,
        inv=inv,
        local_folder=pjoin(save_folder, "event.mseed"),
        host_name=None,
        user_name=None,
        remote_folder=None,
        save_file_name=pjoin(save_folder, "event.mseed"),
    )

    event = events[0]

    """ 2. Read in FIXED parameters """
    if event.baz is not None:
        baz = event.baz
    else:
        baz = input_file["PARAMETERS"]["FIXED"]["baz"]
    if event.distance is not None:
        epi = event.distance
    else:
        epi = input_file["PARAMETERS"]["FIXED"]["epi"]
    if event.origin_time is not None:
        o_time = event.origin_time
    else:
        o_time = input_file["PARAMETERS"]["FIXED"]["o_time"]
    depth = input_file["PARAMETERS"]["FIXED"]["depth"]
    kind = input_file["PARAMETERS"]["FIXED"]["kind"]
    dt = input_file["PARAMETERS"]["FIXED"]["dt"]

    phases = input_file["EVENTS"][event.name]["phases"]
    comps = input_file["EVENTS"][event.name]["components"]
    phase_corrs = input_file["EVENTS"][event.name]["phase_corrs"]
    tstars = input_file["EVENTS"][event.name]["tstars"]
    t_pres = input_file["EVENTS"][event.name]["t_pre"]
    t_posts = input_file["EVENTS"][event.name]["t_post"]
    ylims = input_file["EVENTS"][event.name]["ylims"]

    fmin = input_file["PARAMETERS"]["FILTER"]["fmin"]
    fmax = input_file["PARAMETERS"]["FILTER"]["fmax"]
    zerophase = input_file["PARAMETERS"]["FILTER"]["zerophase"]

    rec_lat = input_file["PARAMETERS"]["RECEIVER"]["latitude"]
    rec_lon = input_file["PARAMETERS"]["RECEIVER"]["longitude"]

    radius = input_file["PARAMETERS"]["PLANET"]["radius"]
    flattening = input_file["PARAMETERS"]["PLANET"]["flattening"]

    """ 3. Read in starting parameters """
    npz_file = input_file["PARAMETERS"]["START_VALUE"]["npz_file"]
    bm_file = input_file["PARAMETERS"]["START_VALUE"]["bm_file"]
    model = TauPyModel(npz_file)
    MOHO = input_file["PARAMETERS"]["START_VALUE"]["MOHO"]
    strike = input_file["PARAMETERS"]["START_VALUE"]["strike"]
    dip = input_file["PARAMETERS"]["START_VALUE"]["dip"]
    rake = input_file["PARAMETERS"]["START_VALUE"]["rake"]
    M0 = input_file["PARAMETERS"]["START_VALUE"]["M0"]

    """ 4. Create dat file """
    """ 4.1 convert strike, dip, rake, M0 into full moment tensor """
    focal_mech0 = np.asarray(_GreensFunctions.convert_SDR(strike, dip, rake, M0)) / 1e14

    """ 4.2 Create the velocity .dat file"""
    bin_path = input_file["BINARY"]["bin_path"]
    # check if folder exists:
    if not exist(save_folder):
        makedirs(save_folder)
    # check if folder is empty
    if not lsdir(save_folder):
        subprocess.call(f"scp {bin_path} .", shell=True, cwd=save_folder)
    bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
    from obspy.geodetics import gps2dist_azimuth

    epi_in_km, _az, _baz = gps2dist_azimuth(
        lat1=event.latitude,
        lon1=event.longitude,
        lat2=rec_lat,
        lon2=rec_lon,
        a=radius,
        f=flattening,
    )

    baz = _baz
    Create_Vmod.create_dat_file(
        src_depth=depth,
        focal_mech=focal_mech0,
        M0=None,
        epi_in_km=epi_in_km,
        baz=baz,
        dt=dt,
        save_path=save_folder,
        bm_file_path=bm_file,
    )

    """ Create observed data array (for the moment all synthetic)"""
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
    st_obs_w, _ = _PreProcess.prepare_event_data(
        event=event,
        phases=phases,
        components=comps,
        slice=True,
        tts=obs_tt,
        t_pre=t_pres,
        t_post=t_posts,
        filter=True,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        noise_level=False,
    )
    st_obs_full, sigmas_noise = _PreProcess.prepare_event_data(
        event=event,
        phases=phases,
        components=comps,
        slice=False,
        filter=True,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        noise_level=False,
    )

    sigma = np.max(sigmas_noise)

    ## Filter
    # st_obs_dt = st_obs_full[0].stats.delta
    # print(st_obs_dt)
    # ff = np.array([0.1, 0.5])
    # st_obs_full = _Gradient.taper_trace(st_obs_full)
    # st_obs_full = _Gradient.bandpass_filt(st_obs_full, ff[0], ff[1], st_obs_dt)

    # st_obs_w = obspy.Stream()
    # for i, phase in enumerate(phases):
    #     tr_full = st_obs_full.select(channel=f"??{comps[i]}")[0].copy()

    #     tr = tr_full.slice(
    #         starttime=event.origin_time + obs_tt[i] - t_pres[i],
    #         endtime=event.origin_time + obs_tt[i] + t_posts[i],
    #     )
    #     st_obs_w += tr
    """ Get the parameter (m) array ready """
    dat_path = pjoin(save_folder, "crfl.dat")
    with open(dat_path, "r+") as f:
        data = f.readlines()
        f.close()
    moment_param = np.array(data[-8].split(), dtype=float)
    moment_param
    """
    4 different cases that you can try:
    1. Changing only the MOHO depth (depth = True, vpvs = False)
    2. Changing the depths only (depth = True, vpvs = False) 
    3. Changing the vpvs only (depth = False, vpvps = True)
    4. Changing the depth and vpvps (depth = True, vpvps = True)
    """
    import pandas as pd

    dat_file = pd.read_csv(
        dat_path,
        skiprows=3,
        skipfooter=11 + 7,
        header=None,
        delim_whitespace=True,
        engine="python",
    )
    # depths = dat_file.values[0:10:2, 0] + 10.0
    # MOHO = 50.0  # depths[5]
    # vp = dat_file.values[0:10:2, 1] + dat_file.values[0:10:2, 1] * 0.01
    # vs = dat_file.values[0:10:2, 3] + dat_file.values[0:10:2, 3] * 0.01
    m0 = np.hstack((moment_param, MOHO))  # Case 1
    # m0 = np.hstack((moment_param, depths)) # Case 2
    # m0 = np.hstack((np.hstack((moment_param, vp)), vs)) # Case 3
    # m0 = np.hstack((np.hstack((np.hstack((moment_param, depths)), vp)), vs)) # Case 4

    depth = True
    vpvs = False

    """ 
    Start the misfit function with your initial model 
    (i.e., it will run the forward model and calculates the misfit)
    """
    n_it = 2
    fac = 1  # factor that you want to multiply with the gradient
    epsilon = 0.01
    xis = np.zeros(n_it)
    dxi_dms = np.zeros((len(m0), n_it))
    m0s = np.zeros((len(m0), n_it + 1))
    m0s[:, 0] = m0

    src_str = _Gradient.SRC_STR(
        binary_file_path=bin_path,
        path_to_dat=save_folder,
        phases=phases,
        components=comps,
        t_pres=t_pres,
        t_posts=t_posts,
        depth=depth,
        dt=dt,
        baz=baz,
        sigmas=sigmas_noise,
        tstars=tstars,
        vpvs=vpvs,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        plot=True,
        st_obs_full=st_obs_full,
        tt_obs=obs_tt,
        ylims=ylims,
    )

    for i in range(n_it):
        xis[i] = src_str.misfit(m0, st_obs_w)
        dxi_dm = af(
            m0,
            src_str.misfit,
            epsilon
            * np.array(
                [
                    np.mean(m0[:-1]),
                    np.mean(m0[:-1]),
                    np.mean(m0[:-1]),
                    np.mean(m0[:-1]),
                    np.mean(m0[:-1]),
                    np.mean(m0[:-1]),
                    m0[-1],
                ]
            ),
            st_obs_w,
        )
        dxi_dms[:, i] = dxi_dm
        m0 -= fac * dxi_dm
        m0s[:, i + 1] = m0
        np.save(pjoin(save_folder, f"it_{i}_dxi_dm.npy"), dxi_dm)
        np.save(pjoin(save_folder, f"it_{i}_xi.npy"), xis[i])
        np.save(pjoin(save_folder, f"it_{i}_m0.npy"), m0)

    np.save(pjoin(save_folder, "dxi_dms.npy"), dxi_dms)
    np.save(pjoin(save_folder, "xis.npy"), xis)
    np.save(pjoin(save_folder, "m0s.npy"), m0s)
