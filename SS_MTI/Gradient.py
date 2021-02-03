"""
Compute gradients of synthetic seismogram w.r.t. model parameters
and returns the misfit directly
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
import scipy.signal as signal
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometer2degrees

import SS_MTI
import Create_Vmod
import subprocess
from SS_MTI import PreProcess as _PreProcess
from SS_MTI import PhaseTracer as _PhaseTracer
from SS_MTI import Misfit as _Misfit
import SS_MTI.SourceTimeFunction as _STF


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


def rot_data(st, angle):
    Z = st.select(channel="xxz")[0].data
    R = st.select(channel="xxr")[0].data
    T = st.select(channel="xxt")[0].data

    data = np.vstack((Z, R, T))

    cost = np.cos(np.deg2rad(angle))
    sint = np.sin(np.deg2rad(angle))
    #     Rot_mat = np.array([[cost, 0, sint], [0, 1, 0], [-sint, 0, cost]]) # ZT (as 2 constant)
    Rot_mat = np.array([[cost, -sint, 0.0], [-sint, cost, 0], [0, 0, 1]])
    rot_data = Rot_mat @ data

    st.select(channel="xxz")[0].data = rot_data[0, :]
    st.select(channel="xxr")[0].data = rot_data[1, :]
    st.select(channel="xxt")[0].data = rot_data[2, :]

    return st


def read_refl_mseeds(path: str, dt: float = None, stack: bool = False):
    """ 
    This functions reads in the output of reflectivity code (miniseeds)
    :param pat:the path where the miniseeds are located
    :param dt: desired dt
    :param stack: stacks the streams in one array (Z,R,T) -> numpy.array
    :returns: obspy stream for stack = False, numpy.array for stack = True
    """
    # path = (
    #     "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/S0235b/jupyter_notebook/235b-location/"
    # )

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

    # st.select(channel="xxz")[0].data *= -1
    # st.select(channel="xxt")[0].data *= -1
    # if st[0].stats.delta != dt:
    #     print(f"Traces will be resampled from dt:{st[0].stats.delta} to dt:{dt}")
    #     st.resample((1 / dt))
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
    ins_full=None,
    ins_w=None,
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
            st_syn_w[i].times() - t_pres[i], -st_syn_w[i].data, lw=3, c="r", label="Reflectivity",
        )
        ax[i].plot(
            st_syn_full[i].times() - tt_syn[i], -st_syn_full[i].data, lw=1, c="r",
        )

        # ax[i].plot(
        #     ins_w[i].times() - t_pres[i], ins_w[i].data, lw=3, c="b", label="Instaseis",
        # )
        # ax[i].plot(
        #     ins_full[i].times() - tt_syn[i], ins_full[i].data, lw=1, c="b",
        # )

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
        ax[i].annotate(
            r"$\times$10$^{%i}$" % (exponent_axis),
            xy=(0.01, 0.75),
            xycoords="axes fraction",
            fontsize=32,
        )
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


def get_tt_from_dat_file(phases: [str], dat_file_path: str, name_tvel: str = None):
    """
    This function can determine travel-times for phases based on crfl.dat file
    (i.e., it extracts the depth, epicentral distance from the .dat file)
    :param phases: phases to determine the traveltimes from
    :param dat_file_path: path to dat_file
    :returns tts: travel times for each phases based on info from .dat file
    """
    if name_tvel is None:
        name_tvel = "temp"

    path_to_tvel = pjoin(dat_file_path, f"{name_tvel}.tvel")
    npz_file = create_taup_model(input_file=path_to_tvel, save_folder=dat_file_path)

    """ 2. Initialize .npz file & calculate phase arrivals """
    Taup = TauPyModel(npz_file)
    depth = Create_Vmod.read_depth_from_dat(dat_file_path)
    epi = Create_Vmod.read_epi_from_dat(dat_file_path)
    tts = []
    for i, phase in enumerate(phases):
        tts.append(_PhaseTracer.get_traveltime(Taup, phase, depth, epi))
    return tts


def window(
    st: obspy.Stream,
    phases: [str],
    components: [str],
    tts: [float],
    t_pres: [float],
    t_posts: [float],
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = False,
):
    """ 
    Windowing an obspy stream around phase arrivals on specific components
    :param st_syn: obspy stream not-windowed
    :param phases: list of phases that you want to window
    :param components: list of components that the phases need to be windowed on 
    :param tts: list of travel-times for each phase
    :returns st_windowed:obspy.Stream, st_full:obspy.Stream, d_stack:np.array 
    """
    assert len(phases) == len(components), "phases and components should be of same length"

    """ Window the synthetic data (& TODO: convolve with STF) """
    d_stack = np.array([])
    st_windowed = obspy.Stream()  # keep this one for plotting
    st_full = obspy.Stream()

    for i, phase in enumerate(phases):
        tr_full = st.select(channel=f"xx{components[i]}")[0].copy()
        """ step 4. convolve with STF """
        # TODO: implement tstars
        # if self.tstars is not None:
        #     stf_len_sec = 30.0
        #     dt_stf = 0.36245500941554226
        #     npts_stf = 82
        #     stf = _STF.stf_tstar(
        #         tstar=self.tstars[i], dt=self.dt, npts=int(stf_len_sec / self.dt)
        #     )[0]
        #     tr_full.data = _STF.convolve_stf(stf, tr_full.data)
        #     stf = _STF.stf_tstar(
        #         tstar=self.tstars[i], dt=db.info.dt, npts=int(stf_len_sec / db.info.dt)
        #     )[0]
        if fmin is not None and fmax is not None:
            _PreProcess.filter_tr(tr_full, fmin=fmin, fmax=fmax, zerophase=zerophase)
        o_time = tr_full.stats.starttime
        tr = tr_full.slice(
            starttime=o_time + tts[i] - t_pres[i], endtime=o_time + tts[i] + t_posts[i],
        )
        st_windowed += tr
        st_full += tr_full
        d_stack = np.hstack((d_stack, tr.data))
    return st_windowed, st_full, d_stack


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
        dt: [float],
        baz: float,
        sigmas: [float],
        tstars: [float] = None,
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
        :param dt: desired delta
        :param baz: backazimuth
        :param sigmas: list of noise level estimate
        :param fmin: lower bound frequency
        :param fmax: upper bound frequency
        :param zerophase: zerophase the data while filtering
        :param plot: plotting or not
        :param st_obs_full: plot = True, observed tream (ordered as phases/component)
        :param tt_obs: observed travel times (ordered as (phases/components)
        :param ylims: limit of the y-axis for plot
        """
        assert len(phases) == len(components), "phases and components should be of same length"
        self.phases = phases
        self.components = components
        self.t_pres = t_pres
        self.t_posts = t_posts
        self.fmin = fmin
        self.fmax = fmax
        self.zerophase = zerophase
        self.dt = dt
        self.tstars = tstars
        self.baz = baz
        self.sigmas = sigmas

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
        return read_refl_mseeds(path=self.f_dat, dt=self.dt, stack=False)

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

        """ Temporary compare with instaseis"""
        # import instaseis
        # db = instaseis.open_db("http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/")
        # # # Event parameters
        # depth_ev = 41.0e3  # depth [m]
        # time_ev = st_syn[0].stats.starttime
        # # Retrieve the waveform
        # receiver = instaseis.Receiver(
        #     latitude=4.5, longitude=135.623447, network="XB", station="ELYSE"
        # )
        # strike = 0
        # dip = 40
        # rake = -150
        # M0 = 9.1925923e13
        # from SS_MTI import GreensFunctions
        # MT = GreensFunctions.convert_SDR(strike, dip, rake, M0)
        # source = instaseis.Source(
        #     latitude=10.99032013,
        #     longitude=160.9467524,
        #     depth_in_m=depth_ev,
        #     m_rr=MT[0],
        #     m_tt=MT[1],
        #     m_pp=MT[2],
        #     m_rt=MT[3],
        #     m_rp=MT[4],
        #     m_tp=MT[5],
        #     origin_time=time_ev,
        #     time_shift=None,
        #     sliprate=None,
        # )
        # st_ins1 = db.get_seismograms(
        #     source=source, receiver=receiver, components="ZRT", dt=0.05, kind="displacement"
        # )

        """ Window the data """
        self.syn_tts = get_tt_from_dat_file(self.phases, self.f_dat, self.tvel_file_name)
        self.st_syn_w, self.st_syn_full, s_syn = window(
            st_syn,
            self.phases,
            self.components,
            self.syn_tts,
            self.t_pres,
            self.t_posts,
            self.fmin,
            self.fmax,
            self.zerophase,
        )

        """ Plot the synthetic data vs. observed data """
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
                ins_w=None,
                ins_full=None,
            )

        """ Get stacked version of observed data """
        s_obs = np.array([])
        for i, phase in enumerate(self.phases):
            s_obs = np.hstack((s_obs, st_obs_w[i].data))

        """ Calculate misfit """
        fmisfit = _Misfit.L2()
        xis = fmisfit.run_misfit(self.phases, st_obs_w, self.st_syn_w, self.sigmas ** 2)
        xi = np.sum(xis)
        # xi = np.linalg.norm(((s_syn - s_obs) / np.max(self.sigmas)), ord=2)
        print(xi)
        self.mi += 1
        if self.mi == len(m) + 2:
            self.mi = 0
            self.it += 1
            if not exist(pjoin(self.f_dat, f"Iteration_{self.it}")):
                makedirs(pjoin(self.f_dat, f"Iteration_{self.it}"))
        return xi
