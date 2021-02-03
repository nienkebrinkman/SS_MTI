"""
Interactive example to invert focal mechanism and structure simultaneously
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
from obspy.signal.cross_correlation import xcorr_max, correlate
import scipy.signal as signal
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from scipy.optimize import approx_fprime as af
import re

import SS_MTI
import Create_Vmod
from EventInterface import EventObj
import subprocess
from SS_MTI import PostProcessing as _PostProcessing
from SS_MTI import PreProcess as _PreProcess
from SS_MTI import PhaseTracer as _PhaseTracer
from SS_MTI import Gradient as _Gradient

# save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/m0_gradient_descent"
# grad = np.load(pjoin(save_path, "it_1_dxi_dm.npy"))
# misfit = np.load(pjoin(save_path, "it_1_xi.npy"))
# m0s = np.load(pjoin(save_path, "it_0_m0.npy"))
# m1s = np.load(pjoin(save_path, "it_1_m0.npy"))

# a = 1


def create_taup_model(input_file: str, save_folder: str, ret=True):
    """
    creates .npz file that can be used for obspy ray-tracing
    :param input_file: Either .nd or .tvel file
    :param save_folder: folder where the .npz file will be saved
    :param ret: if True, return value will be given (False = No returns)
    :returns full .npz file path
    """
    from obspy.taup.taup_create import build_taup_model

    build_taup_model(input_file, save_folder)

    if ret:
        return input_file.replace(".tvel", ".npz")


def Get_location(la_s, lo_s, la_r, lo_r, radius=3389.5, flattening=0):
    """
    Get the epicentral distance, azimuth and backazimuth
    """
    dist, az, baz = gps2dist_azimuth(
        lat1=la_s, lon1=lo_s, lat2=la_r, lon2=lo_r, a=radius, f=flattening
    )
    epi = kilometer2degrees(dist, radius=radius)
    return epi, az, baz


# path of the reclectivity binary:
bin_path = (
    "/home/nienke/Documents/Research/SS_MTI/External_packages/reflectivity_Mars/SRC/test/crfl_sac"
)

""" Define starting parameters """
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 0
lon_src = 0
depth = 20.0
name = "Test_Event"

npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
model = TauPyModel(npz_file)

lat_rec = -30
lon_rec = 0

dt = 0.025

epi, az, baz = Get_location(lat_src, lon_src, lat_rec, lon_rec)
print(epi)
epi_in_km = 1774.7380

m_rr0 = 0.2
m_tt0 = 0.8
m_pp0 = 0.0
m_rt0 = 0.0
m_rp0 = 0.0
m_tp0 = 0.0
focal_mech0 = [m_rr0, m_tt0, m_pp0, m_rt0, m_rp0, m_tp0]


# This is basically your prior model, you need to set this up once:
save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/m0_gradient_descent/"


""" Create observed data array (for the moment all synthetic)"""
path_observed = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/obs/"
st_obs = _Gradient.read_refl_mseeds(path_observed, stack=False)
""" Travel times of observed data """
npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
Taup = TauPyModel(npz_file)
obs_tts = []
for i, phase in enumerate(phases):
    obs_tts.append(_PhaseTracer.get_traveltime(Taup, phase, depth, epi))
""" Window the observed data """
st_obs_w = obspy.Stream()  # keep this one for plotting
st_obs_full = obspy.Stream()
for i, phase in enumerate(phases):
    tr_full = st_obs.select(channel=f"xx{comps[i]}")[0].copy()
    if fmin is not None and fmax is not None:
        _PreProcess.filter_tr(tr_full, fmin=fmin, fmax=fmax, zerophase=zerophase)
    o_time = tr_full.stats.starttime
    tr = tr_full.slice(
        starttime=o_time + obs_tts[i] - t_pres[i], endtime=o_time + obs_tts[i] + t_posts[i],
    )
    st_obs_w += tr
    st_obs_full += tr_full


""" """
# check if folder exists:
if not exist(save_path):
    makedirs(save_path)
# check if folder is empty
if not lsdir(save_path):
    subprocess.call(f"scp {bin_path} .", shell=True, cwd=save_path)
bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
Create_Vmod.create_dat_file(
    src_depth=depth,
    focal_mech=focal_mech0,
    M0=None,
    epi_in_km=epi_in_km,
    baz=baz,
    dt=dt,
    save_path=save_path,
    bm_file_path=bm_file_path,
)

""" Get the parameter (m) array ready """
import pandas as pd

dat_path = pjoin(save_path, "crfl.dat")
dat_file = pd.read_csv(
    dat_path, skiprows=3, skipfooter=11 + 7, header=None, delim_whitespace=True, engine="python",
)

with open(dat_path, "r+") as f:
    data = f.readlines()
    f.close()
moment_param = np.array(data[-8].split(), dtype=float)
# depths = dat_file.values[0:10:2, 0] + 10.0
MOHO = 77.3680
vp = dat_file.values[0:10:2, 1] + dat_file.values[0:10:2, 1] * 0.01
vs = dat_file.values[0:10:2, 3] + dat_file.values[0:10:2, 3] * 0.01

"""
4 different cases that you can try:
1. Changing only the MOHO depth (depth = True, vpvs = False)
2. Changing the depths only (depth = True, vpvs = False) 
3. Changing the vpvs only (depth = False, vpvps = True)
4. Changing the depth and vpvps (depth = True, vpvps = True)
"""
m0 = np.hstack((moment_param, MOHO))  # Case 1
# m0 = np.hstack((moment_param, depths)) # Case 2
# m0 = np.hstack((np.hstack((moment_param, vp)), vs)) # Case 3
# m0 = np.hstack((np.hstack((np.hstack((moment_param, depths)), vp)), vs)) # Case 4

depth = True
vpvs = False

sigmas = np.ones(len(phases)) * 1e-10

""" 
Start the misfit function with your initial model 
(i.e., it will run the forward model and calculates the misfit)
"""
n_it = 15
fac = 0.001  # factor that you want to multiply with the gradient
epsilon = 0.01
xis = np.zeros(n_it)
dxi_dms = np.zeros((len(m0), n_it))
m0s = np.zeros((len(m0), n_it + 1))
m0s[:, 0] = m0

src_str = _Gradient.SRC_STR(
    binary_file_path=bin_path,
    path_to_dat=save_path,
    phases=phases,
    components=comps,
    t_pres=t_pres,
    t_posts=t_posts,
    depth=depth,
    vpvs=vpvs,
    fmin=fmin,
    fmax=fmax,
    dt=dt,
    baz=baz,
    sigmas=sigmas,
    zerophase=zerophase,
    plot=True,
    st_obs_full=st_obs_full,
    tt_obs=obs_tts,
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
    np.save(pjoin(save_path, "dxi_dms.npy"), dxi_dms)
    np.save(pjoin(save_path, "xis.npy"), xis)
    np.save(pjoin(save_path, "m0s.npy"), m0s)

# np.save(pjoin(save_path, "dxi_dms.npy"), dxi_dms)
# np.save(pjoin(save_path, "xis.npy"), xis)
# np.save(pjoin(save_path, "m0s.npy"), m0s)
