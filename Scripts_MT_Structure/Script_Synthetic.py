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

import SS_MTI
import Create_Vmod
from EventInterface import EventObj
import subprocess
from SS_MTI import PostProcessing as _PostProcessing


def Get_location(la_s, lo_s, la_r, lo_r, radius=3389.5, flattening=0):
    """
    Get the epicentral distance, azimuth and backazimuth
    """
    dist, az, baz = gps2dist_azimuth(
        lat1=la_s, lon1=lo_s, lat2=la_r, lon2=lo_r, a=radius, f=flattening
    )
    epi = kilometer2degrees(dist, radius=radius)
    return epi, az, baz


def bandpass_filt(stream, f1, f2, dt, order=2):
    """
    filter waveforms by using Butterworth bandpass filter
    """
    fN = 1.0 / (2.0 * dt)
    low = f1 / fN
    high = f2 / fN
    [b, a] = signal.butter(order, low, "high")
    [d, c] = signal.butter(order, high, "low")

    for tr in stream:
        tr.data = signal.filtfilt(b, a, tr.data)
        tr.data = signal.filtfilt(d, c, tr.data)

    return stream


def taper_trace(stream, max_perc=0.03, ttype="hann"):
    """
    taper traces in stream
    """
    for tr in stream:
        tr.taper(max_perc, ttype)

    return stream


def scale_traces(st1, st2):
    """
    scale each trace in st2 to the maximum value of each trace in st1
    """
    for i in range(len(st1)):
        max_val = np.max(st1[i].data)
        if np.max(st2[i].data) != 0:
            st2[i].data = st2[i].data / np.max(st2[i].data) * max_val

    return st1, st2


def read_refl_mseeds(path: str, stack: bool = False):
    """ 
    This functions reads in the output of reflectivity code (miniseeds)
    :param pat:the path where the miniseeds are located
    :param stack: stacks the streams in one array (Z,R,T) -> numpy.array
    :returns: obspy stream for stack = False, numpy.array for stack = True
    """

    st = obspy.Stream()
    st_files = [f for f in lsdir(path) if f.startswith("st") if isfile(pjoin(path, f))]

    for st_file in st_files:
        st_temp = obspy.read(pjoin(path, st_file))

        st_temp[0].stats.channel = "xx" + st_file.split(".")[-1]
        distance = st_temp[0].stats.sac.gcarc
        st_temp[0].stats.distance = distance

        B = st_temp[0].stats.sac.b  # beginning time
        tstart = st_temp[0].stats.starttime  # absolute starttime
        et0 = tstart - B
        st_temp[0] = st_temp[0].trim(
            starttime=tstart - B, endtime=st_temp[0].stats.endtime, pad=True, fill_value=0.0
        )
        st += st_temp

    st.select(channel="xxz")[0].data *= -1
    st.select(channel="xxt")[0].data *= -1
    if stack:
        Z = st.select(channel="xxz")[0].data
        R = st.select(channel="xxr")[0].data
        T = st.select(channel="xxt")[0].data
        return np.hstack((np.hstack((Z, R)), T))
    else:
        return st


class SRC_STR:
    """ 
    This class does source and structure inversion
    """

    def __init__(self, binary_file_path: str, path_to_dat: str, vpvs: bool, depth: bool):
        """ 
        If vpvs and depth are both True (depth and vpvs are both inverted)
        :param binary_path: path to reflectivity code binary file
        :param path_to_dat: path to dat file
        :param vpvs: if true vpvs updates will be done
        :param depth: if true layer depth updates will be done 
        """

        """ Step 1. copy binary file into datfile folder"""
        subprocess.call(f"scp {binary_file_path} .", shell=True, cwd=path_to_dat)

        self.f_dat = path_to_dat
        self.it = 0
        self.vpvs = vpvs
        self.depth = depth

    def forward(self, m: np.array):
        """
        forward function that runs the reflectivity code based on the model parameters(m).
        :param m: array of source and structure parameters
        :returns s_syn: synthetic output based on input m
        """

        """ step 1. plug in the model parameters in the .dat file """
        Create_Vmod.update_dat_file(self.f_dat, m, self.vpvs, self.depth)
        """ step 2. run the dat file """
        # subprocess.call("./crfl_sac", shell=True, cwd=self.f_dat)
        """ step 3. read in the output of the run """
        return read_refl_mseeds(path=self.f_dat, stack=True)

    def misfit(self, m: np.array, s_obs):
        """
        Misfit function (L2) based on model parameters (m) and observed data.
        It runs the forward model (reflectvitiy code under the hood)
        :param m: array of source and structure parameters
        :param s_obs: array with observed data
        :returns xi: misfit of synthetic data vs observed data
        """
        print(f"model paramer nr {self.it} is updated now")
        """ Run the forward modeller"""
        s_syn = self.forward(m)

        """ Windowing s_syn and s_obs around P and S"""

        plt.plot(s_obs, c="r", label="obs")
        plt.plot(s_syn, c="b", label="syn")
        plt.legend()
        plt.savefig(pjoin(self.f_dat, f"syn_obs_{self.it}.svg"))
        plt.close()

        xi = np.linalg.norm((s_syn - s_obs), ord=2)
        self.it += 1
        return xi


# path of the reclectivity binary:
bin_path = (
    "/home/nienke/Documents/Research/SS_MTI/External_packages/reflectivity_Mars/SRC/test/crfl_sac"
)

""" Create observed data array """
path_observed = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/obs/"
s_obs = read_refl_mseeds(path_observed, stack=True)

""" Define starting parameters """
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 0
lon_src = 0
depth = 20.0
name = "Test_Event"

npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
model = TauPyModel(npz_file)

lat_rec = -20
lon_rec = 0

dt = 0.025

epi, az, baz = Get_location(lat_src, lon_src, lat_rec, lon_rec)

m_rr0 = 0.0
m_tt0 = 1.0
m_pp0 = 0.0
m_rt0 = 0.0
m_rp0 = 0.0
m_tp0 = 0.0
focal_mech0 = [m_rr0, m_tt0, m_pp0, m_rt0, m_rp0, m_tp0]

# This is basically your prior model, you need to set this up once:
save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/m0/"
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
    epi=epi,
    baz=baz,
    save_path=save_path,
    bm_file_path=bm_file_path,
)

""" Get the parameter (m) array ready """
import pandas as pd

dat_path = pjoin(save_path, "crfl.dat")
dat_file = pd.read_csv(
    dat_path, skiprows=3, skipfooter=11 + 7, header=None, delim_whitespace=True, engine="python",
)
moment_param = np.loadtxt(dat_path, skiprows=103, max_rows=1)
depths = dat_file.values[0:10:2, 0] + 10.0
MOHO = 40.0  # depths[5]
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

""" 
Start the misfit function with your initial model 
(i.e., it will run the forward model and calculates the misfit)
"""
n_it = 1
fac = 1  # factor that you want to multiply with the gradient
xis = np.zeros(n_it)
dxi_dms = np.zeros((len(m0), n_it))
m0s = np.zeros((len(m0), n_it + 1))
m0s[:, 0] = m0
for i in range(n_it):
    src_str = SRC_STR(binary_file_path=bin_path, path_to_dat=save_path, depth=depth, vpvs=vpvs)
    xis[i] = src_str.misfit(m0, s_obs)
    dxi_dm = af(m0, src_str.misfit, 0.01, s_obs)
    dxi_dms[:, i] = dxi_dm
    m0 -= fac * dxi_dm
    m0s[:, i + 1] = m0

np.save(pjoin(save_path, "dxi_dms.npy"), dxi_dms)
np.save(pjoin(save_path, "xis.npy"), xis)
np.save(pjoin(save_path, "m0s.npy"), m0s)
