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


def plot_waveforms_src_str(
    save_folder: str,
    st_syn_full: obspy.Stream,
    st_obs_full: obspy.Stream,
    st_syn_w: obspy.Stream,
    st_obs_w: obspy.Stream,
    phases: [str],
    components: [str],
    t_pres: [float],
    t_posts: [float],
    tt_obs: [float],
    tt_syn: [float],
    ylims: [float],
    save_file: str,
):
    """
    plot function

    :param st_syn_full: synthetic obspy stream (not windowed)
    :param st_obs_full: observed obspy stream (not windowed)
    :param st_syn_w: synthetic obspy stream (windowed)
    :param st_obs_w: observed obspy stream (windowed)
    :param phases: phases that are windowed
    :param components: on which component the phases should be windowed
    :param t_pres: time before phase arrival to window
    :param t_posts: time after phase arrival to window
    :param tt_obs: traveltime arrivals of observed data
    :param tt_syn: traveltime arrivals of synthetic data
    :param ylims: ylimits (ordered as in phases)
    :param save_file: path and name of the file 
    """
    fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(18, 20))
    for i, phase in enumerate(phases):
        ax[i].plot(
            st_obs_w[i].times() - t_pres[i], st_obs_w[i].data, lw=3, c="k", label="Observed",
        )
        ax[i].plot(
            st_obs_full[i].times() - tt_obs[i], st_obs_full[i].data, lw=1, c="k",
        )

        ax[i].plot(
            st_syn_w[i].times() - t_pres[i], st_syn_w[i].data, lw=3, c="r", label="Synthetic",
        )
        ax[i].plot(
            st_syn_full[i].times() - tt_syn[i], st_syn_full[i].data, lw=1, c="r",
        )

        ax[i].text(
            s=f"{phase}{components[i]}",
            x=0.99,
            y=0.75,
            ha="right",
            transform=ax[i].transAxes,
            color="blue",
            fontsize=40,
        )
        ax[i].tick_params(axis="both", which="major", labelsize=35)
        ax[i].tick_params(axis="both", which="minor", labelsize=25)
        ax[i].get_yaxis().get_offset_text().set_visible(False)
        ax_max = max(ax[i].get_yticks())
        exponent_axis = np.floor(np.log10(ax_max)).astype(int)
        global_max = max([tr.data.max() for tr in st_obs_w]) * 1.2
        global_min = min([tr.data.min() for tr in st_obs_w]) * 1.2
        if ylims is None:
            ax[i].set_ylim(global_min, global_max)
        else:
            ax[i].set_ylim(-ylims[i], ylims[i])
        ymax = ax[i].get_ylim()[1]
        ax[i].axvline(x=t_posts[i], c="grey", ls="dashed")
        ax[i].axvline(x=-t_pres[i], c="grey", ls="dashed")
        ax[i].axvline(x=0.0, c="dimgrey", lw=2)
        ax[i].text(
            0 + 0.1, ymax * 0.8, phase, verticalalignment="center", color="dimgray", fontsize=30,
        )

    fig.text(0.01, 0.5, "Displacement (nm)", va="center", rotation="vertical", fontsize=45)
    fig.text(
        0.5, 0.88, "Test", ha="center", va="bottom", size="x-large", color="blue", fontsize=45,
    )

    ax[0].legend(
        prop={"size": 35},
        loc="center left",
        bbox_to_anchor=(0.12, 0.93),
        bbox_transform=fig.transFigure,
    )

    ax[-1].set_xlim(-10.0, 32.0)
    ax[-1].set_xlabel("time after phase (s)", fontsize=45)
    plt.savefig(save_file)


class SRC_STR:
    """ 
    This class does source and structure inversion
    """

    def __init__(
        self,
        binary_file_path: str,
        path_to_dat: str,
        phases: [str],
        components: [str],
        t_pres: [float],
        t_posts: [float],
        vpvs: bool,
        depth: bool,
        fmin: float = None,
        fmax: float = None,
        zerophase: bool = False,
        plot: bool = True,
        st_obs_full: obspy.Stream = None,
        tt_obs: [float] = None,
        ylims: [float] = None,
    ):
        """ 
        If vpvs and depth are both True (depth and vpvs are both inverted)
        :param binary_path: path to reflectivity code binary file
        :param path_to_dat: path to dat file
        :param phases: phases to be windowed
        :param components: on which component the phases should be windowed
        :param t_pres: time before phase arrival to window
        :param t_posts: time after phase arrival to window
        :param vpvs: if true vpvs updates will be done
        :param depth: if true layer depth updates will be done 
        :param fmin: lower bound frequency
        :param fmax: upper bound frequency
        :param zerophase: zerophase the data while filtering
        :param plot: plotting or not
        :param st_obs_full: plot = True, observed tream (ordered as phases/component)
        :param tt_obs: observed travel times (ordered as (phases/components)
        """
        assert len(phases) == len(components), "phases and components should be of same length"
        self.phases = phases
        self.components = components
        self.t_pres = t_pres
        self.t_posts = t_posts
        self.fmin = fmin
        self.fmax = fmax
        self.zerophase = zerophase

        """ Step 1. copy binary file into datfile folder"""
        subprocess.call(f"scp {binary_file_path} .", shell=True, cwd=path_to_dat)

        self.f_dat = path_to_dat
        self.mi = 0  # model parameter that is currently updated
        self.it = 0  # iteration of the full inversion
        if not exist(pjoin(self.f_dat, f"Iteration_{self.it}")):
            makedirs(pjoin(self.f_dat, f"Iteration_{self.it}"))
        self.vpvs = vpvs
        self.depth = depth

        """ Plotting values """
        self.plot = plot
        if self.plot:
            assert (
                st_obs_full is not None and tt_obs is not None
            ), "if plot = True, you must give st_obs_full"
            self.st_obs_full = st_obs_full
            self.tt_obs = tt_obs
            self.ylims = ylims

    def forward(self, m: np.array):
        """
        forward function that runs the reflectivity code based on the model parameters(m).
        :param m: array of source and structure parameters
        :returns s_syn: synthetic output based on input m
        """

        """ step 1. plug in the model parameters in the .dat file (.tvel will be created)"""
        create_tvel = True
        self.tvel_file_name = f"{self.mi}"
        Create_Vmod.update_dat_file(
            self.f_dat, m, self.vpvs, self.depth, create_tvel, self.tvel_file_name
        )
        """ step 2. run the dat file """
        subprocess.call("./crfl_sac", shell=True, cwd=self.f_dat)
        """ step 3. read in the output of the run """
        return read_refl_mseeds(path=self.f_dat, stack=False)

    def misfit(self, m: np.array, st_obs_w: obspy.Stream):
        """
        Misfit function (L2) based on model parameters (m) and observed data.
        It runs the forward model (reflectvitiy code under the hood)
        :param m: array of source and structure parameters
        :param st_obs_w: obspy stream (windowed) with observed data
        :returns xi: misfit of synthetic data vs observed data
        """
        print(f"model parameter nr {self.mi} is updated now")
        """ Run the forward modeller"""
        st_syn = self.forward(m)

        """ Windowing s_syn and s_obs around P and S """
        """ 1. .npz file needs to be created from .tvel file """
        path_to_tvel = pjoin(self.f_dat, self.tvel_file_name + ".tvel")
        npz_file = create_taup_model(input_file=path_to_tvel, save_folder=self.f_dat)

        """ 2. Initialize .npz file & calculate phase arrivals """
        Taup = TauPyModel(npz_file)
        depth = Create_Vmod.read_depth_from_dat(self.f_dat)
        epi = Create_Vmod.read_epi_from_dat(self.f_dat)
        self.syn_tts = []
        for i, phase in enumerate(self.phases):
            self.syn_tts.append(_PhaseTracer.get_traveltime(Taup, phase, depth, epi))

        """ 3. Window the synthetic data """
        s_syn = np.array([])
        s_obs = np.array([])
        self.st_syn_w = obspy.Stream()  # keep this one for plotting
        self.st_syn_full = obspy.Stream()
        for i, phase in enumerate(self.phases):
            tr_full = st_syn.select(channel=f"xx{self.components[i]}")[0].copy()
            if self.fmin is not None and self.fmax is not None:
                _PreProcess.filter_tr(
                    tr_full, fmin=self.fmin, fmax=self.fmax, zerophase=self.zerophase
                )
            o_time = tr_full.stats.starttime
            tr = tr_full.slice(
                starttime=o_time + self.syn_tts[i] - self.t_pres[i],
                endtime=o_time + self.syn_tts[i] + self.t_posts[i],
            )
            self.st_syn_w += tr
            self.st_syn_full += tr_full
            s_syn = np.hstack((s_syn, tr.data))
            s_obs = np.hstack((s_obs, st_obs_w.select(channel=f"xx{self.components[i]}")[0].data))

        """ 4. Plot the synthetic data vs. observed data """
        if self.plot:
            plot_waveforms_src_str(
                save_folder=self.f_dat,
                st_syn_full=self.st_syn_full,
                st_obs_full=self.st_obs_full,
                st_syn_w=self.st_syn_w,
                st_obs_w=st_obs_w,
                phases=self.phases,
                components=self.components,
                t_pres=self.t_pres,
                t_posts=self.t_posts,
                tt_obs=self.tt_obs,
                tt_syn=self.syn_tts,
                ylims=self.ylims,
                save_file=pjoin(self.f_dat, f"Iteration_{self.it}", f"Update_{self.mi}.svg"),
            )

        """ 5. Calculate misfit """
        xi = np.linalg.norm((s_syn - s_obs), ord=2)
        print(xi)
        self.mi += 1
        if self.mi == len(m0) + 2:
            self.mi = 0
            self.it += 1
            if not exist(pjoin(self.f_dat, f"Iteration_{self.it}")):
                makedirs(pjoin(self.f_dat, f"Iteration_{self.it}"))
        return xi


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

m_rr0 = 0.0
m_tt0 = 1.0
m_pp0 = 0.0
m_rt0 = 0.0
m_rp0 = 0.0
m_tp0 = 0.0
focal_mech0 = [m_rr0, m_tt0, m_pp0, m_rt0, m_rp0, m_tp0]

phases = ["P", "P", "S", "S", "S"]
comps = ["Z", "R", "Z", "R", "T"]
t_pres = [1, 1, 1, 1, 1]
t_posts = [30, 30, 30, 30, 30]
ylims = [1e-9, 1e-9, 2e-9, 3e-9, 2e-9]

fmin = 0.2
fmax = 0.6
zerophase = False

# This is basically your prior model, you need to set this up once:
save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/m0/"


""" Create observed data array (for the moment all synthetic)"""
path_observed = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/obs/"
st_obs = read_refl_mseeds(path_observed, stack=False)
""" Travel times of observed data """
npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz"
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

with open(dat_path, "r+") as f:
    data = f.readlines()
    f.close()
moment_param = np.array(re.findall("\d+\.\d+", data[-8]), dtype=float)
depths = dat_file.values[0:10:2, 0] + 10.0
MOHO = 50.0  # depths[5]
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
    src_str = SRC_STR(
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
        zerophase=zerophase,
        plot=True,
        st_obs_full=st_obs_full,
        tt_obs=obs_tts,
        ylims=ylims,
    )

    xis[i] = src_str.misfit(m0, st_obs_w)
    dxi_dm = af(m0, src_str.misfit, 0.01, st_obs_w)
    dxi_dms[:, i] = dxi_dm
    m0 -= fac * dxi_dm
    m0s[:, i + 1] = m0

np.save(pjoin(save_path, "dxi_dms.npy"), dxi_dms)
np.save(pjoin(save_path, "xis.npy"), xis)
np.save(pjoin(save_path, "m0s.npy"), m0s)
