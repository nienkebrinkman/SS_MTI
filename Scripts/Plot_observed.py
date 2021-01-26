"""
Plot observed data (time domain) and spectra of pre-event and event
Useful to determine the bandpass filter
Plots 3:
1. Original data in time and frequency domain
2. P and S wave windows in time domain
3. Hodograms of P and S waves
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from os import listdir as lsdir
import instaseis
import numpy as np
import matplotlib.pyplot as plt
import obspy
from obspy import UTCDateTime as utct

import SS_MTI
import EventInterface
from SS_MTI import PostProcessing as _PostProcessing
from SS_MTI import PreProcess as _PreProcess

# Event325a:
# t_start = -700
# t_end = 1100
# pre_noise = 600.0  # Time for the noise window  (for event S0409d: 1100)
# S_zoom = 5e-9
# P_zoom = 1e-9
# P_pre = 5.0
# P_post = 5.0
# S_pre = 10.0
# S_post = 10.0
# P_shift = 0.0
# S_shift = 5.0
# left_BP = 0.2
# right_BP = 0.6

# Event 407a:
# t_start = -700
# t_end = 1100
# pre_noise = 600.0  # Time for the noise window  (for event S0409d: 1100)
# S_zoom = 2e-10
# P_zoom = 2e-10
# P_shift = 0.0
# S_shift = 1.7
# P_pre = 5.0
# P_post = 5.0
# S_pre = 5.0
# S_post = 5.0
# left_BP = 0.3
# right_BP = 0.7

Event_names = {"S0409d": {}}
components = ["Z", "N", "E"]
save_folder = "/home/nienke/Documents/Research/Data/MTI/Data_2021/"

""" Parameters for Time domain plot (set these parameters to None if you dont know) """
t_start = -1200
t_end = 1100
pre_noise = 1100.0  # Time for the noise window  (for event S0409d: 1100)
S_zoom = 2e-10
P_zoom = 2e-10

P_shift = 0.0
S_shift = 0.0
""" This is for the Hodogram plot: """
P_pre = 1.0
P_post = 20.0
S_pre = 5.0
S_post = 5.0

""" Parameters for Frequency domain plot (set these parameters to None if you dont know) """
left_BP = 0.2
right_BP = 0.8


""" Read the inventory and catalog file (the once that contain info about the marsquakes) """
# path = "/home/nienke/Documents/Research/Data/MTI/Catalog/old_catalog"
path = "/home/nienke/Documents/Research/Data/MTI/Catalog/"
# path = "/home/nienke/Documents/Research/SS_MTI/Data/"
path_to_inventory = pjoin(path, "inventory.xml")
path_to_catalog = pjoin(path, "catalog.xml")
inv = SS_MTI.DataGetter.read_inv(inv_path=path_to_inventory)  # Inventory file
cat = SS_MTI.DataGetter.read_cat(cat_path=path_to_catalog)  # Catalog file
# for ev in cat:
#     print(ev.name)

""" Collect the corresponding data that belongs to the events """
events = SS_MTI.DataGetter.read_events_from_cat(
    event_params=Event_names,
    cat=cat,
    inv=inv,
    local_folder="/mnt/marshost/",
    host_name="marshost.ethz.ch",
    user_name="sysop",
    remote_folder="/data/",
    save_file_name=pjoin(save_folder, "event.mseed"),
)

""" Specify receiver """
lat_rec = 4.5  # 02384
lon_rec = 135.623447
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

for event in events:
    st = event.waveforms_VBB.copy()
    # epi, az, baz = _PreProcess.Get_location(
    #     la_s=event.latitude, lo_s=event.longitude, la_r=rec.latitude, lo_r=rec.longitude
    # )
    # st = st.rotate("NE->RT", back_azimuth=baz)

    P_arr = utct(event.picks["P"]) - utct(event.origin_time) + P_shift
    S_arr = utct(event.picks["S"]) - utct(event.origin_time) + S_shift

    """ Plot the waveforms and spectra """
    fig, ax = plt.subplots(nrows=3, ncols=2, sharex="col", sharey="col", figsize=(15, 15))
    for i, comp in enumerate(components):
        tr = st[i].copy()
        o_time = event.origin_time - tr.stats.starttime

        """ Plot the origin time, P and S arrivals """
        ax[i, 0].axvline(x=0, c="gold", ls="dashed", label="Origin time")
        ax[i, 0].axvline(x=P_arr, c="red", ls="dashed", label="P-arrival")
        ax[i, 0].axvline(x=S_arr, c="blue", ls="dashed", label="S-arrival")

        """ Plot the original data in time-domain """
        ax[i, 0].plot(
            tr.times() - o_time, tr.data, label=f"Raw data", c="k",
        )

        """ Slice data """
        tr_event = tr.slice(
            starttime=event.origin_time + P_arr - 5.0, endtime=event.origin_time + S_arr + 30.0,
        )
        tr_P = tr.slice(
            starttime=event.origin_time + P_arr - 5.0, endtime=event.origin_time + P_arr + 50.0,
        )
        tr_S = tr.slice(
            starttime=event.origin_time + S_arr - 5.0, endtime=event.origin_time + S_arr + 70.0,
        )
        tr_noise = tr.slice(
            starttime=event.origin_time - pre_noise, endtime=event.origin_time - pre_noise + 200.0,
        )

        """ Plot sliced data """
        ax[i, 0].plot(
            tr_event.times() + P_arr - 5.0, tr_event.data, label=f"selected event data", c="gold",
        )
        ax[i, 0].plot(
            tr_P.times() + P_arr - 5.0, tr_P.data, label=f"selected P data", c="red",
        )
        ax[i, 0].plot(
            tr_S.times() + S_arr - 5.0, tr_S.data, label=f"selected S data", c="blue",
        )
        ax[i, 0].plot(
            tr_noise.times() - pre_noise,
            tr_noise.data,
            label=f"selected noise data",
            c="darkgreen",
        )

        """ Convert data to freq-domain """
        f_event, p_event = _PreProcess.calc_PSD(tr_event, winlen_sec=20.0)
        f_P, p_P = _PreProcess.calc_PSD(tr_P, winlen_sec=20.0)
        f_S, p_S = _PreProcess.calc_PSD(tr_S, winlen_sec=20.0)
        f_noise, p_noise = _PreProcess.calc_PSD(tr_noise, winlen_sec=20.0)

        """ Plot various t-star values """
        colors = ["violet", "darkred", "limegreen", "pink"]
        for j, t in enumerate([0.1, 0.5, 1.0, 2.0]):
            p_tstar = p_P[0] * np.exp(-np.pi * f_event * t * 2)
            ax[i, 1].semilogx(
                f_event, 10 * np.log10(p_tstar), label=f"t*:{t}", lw=3, alpha=0.2, c=colors[j]
            )

        """ Plot the original data in Frequency-domain """
        ax[i, 1].semilogx(
            f_event, 10 * np.log10(p_event), c="gold",
        )
        ax[i, 1].semilogx(
            f_P, 10 * np.log10(p_P), c="red",
        )
        ax[i, 1].semilogx(
            f_S, 10 * np.log10(p_S), c="blue",
        )
        ax[i, 1].semilogx(
            f_noise, 10 * np.log10(p_noise), c="darkgreen",
        )

        """ Plot the possible bandpass filter """
        if left_BP is not None and right_BP is not None:
            ax[i, 1].axvspan(left_BP, right_BP, facecolor="gray", alpha=0.2)

        """ Set plotting parameters """
        fsize = 15
        ax[i, 0].text(
            s=comp,
            x=0.99,
            y=0.90,
            ha="right",
            transform=ax[i, 0].transAxes,
            color="blue",
            fontsize=fsize,
        )
        ax[i, 0].set_ylabel("Displacement (m)", fontsize=fsize)
        ax[i, 1].set_ylabel("Displacement PSD (dB)", fontsize=fsize)
        # ax[i, 0].set_xlim(P_arr - 150.0, S_arr + 300.0)
        ax[0, 1].legend(fontsize=fsize, ncol=1)
        ax[0, 0].legend(
            ncol=2,
            prop={"size": fsize},
            loc="center left",
            bbox_to_anchor=(0.6, 0.94),
            bbox_transform=fig.transFigure,
        )
        ax[i, 0].set_ylim(-0.2e-7, 0.2e-7)
        if t_start is not None and t_end is not None:
            ax[i, 0].set_xlim(t_start, t_end)
        ax[i, 1].set_xlim(9e-2, 4e0)
        ax[i, 1].set_ylim(-240, -170)

        ax[-1, 0].set_xlabel("Time (s)", fontsize=fsize)
        ax[-1, 1].set_xlabel("Frequency (Hz)", fontsize=fsize)
    fig.text(
        0.55, 0.9, event.name, ha="right", va="bottom", size="medium", color="black", fontsize=30,
    )
    plt.savefig(pjoin(save_folder, f"{event.name}.svg"))

    """ Plot P and S separate for each component """
    fig, ax = plt.subplots(nrows=3, ncols=2, sharex="col", sharey="col", figsize=(30, 15))

    for i, comp in enumerate(components):
        tr = st[i].copy()
        if left_BP is not None and right_BP is not None:
            _PreProcess.filter_tr(tr, fmin=left_BP, fmax=right_BP, zerophase=False)
        o_time = event.origin_time - tr.stats.starttime

        """ Plot the origin time, P and S arrivals """

        ls = "dashed"
        if P_shift != 0.0:
            ax[i, 0].axvline(x=-P_shift, c="black", ls=ls, label="P-arrival (original)")
            ls = "dotted"
        ax[i, 0].axvline(x=0, c="black", ls=ls, label="P-arrival")
        ls = "dashed"
        if S_shift != 0.0:
            ax[i, 1].axvline(x=-S_shift, c="black", ls=ls, label="S-arrival (original)")
            ls = "dotted"
        ax[i, 1].axvline(x=0, c="black", ls=ls, label="S-arrival")

        """ Slice data """
        tr_P = tr.slice(
            starttime=event.origin_time + P_arr - 20.0, endtime=event.origin_time + P_arr + 50.0,
        )
        tr_S = tr.slice(
            starttime=event.origin_time + S_arr - 20.0, endtime=event.origin_time + S_arr + 70.0,
        )

        """ Plot sliced data """
        ax[i, 0].plot(
            tr_P.times() - 20.0, tr_P.data, label=f"selected P data", c="red",
        )
        ax[i, 1].plot(
            tr_S.times() - 20.0, tr_S.data, label=f"selected S data", c="blue",
        )

        """ Set plotting parameters """
        fsize = 15
        ax[i, 0].text(
            s=f"P{comp}",
            x=0.99,
            y=0.90,
            ha="right",
            transform=ax[i, 0].transAxes,
            color="blue",
            fontsize=fsize,
        )
        ax[i, 1].text(
            s=f"S{comp}",
            x=0.99,
            y=0.90,
            ha="right",
            transform=ax[i, 1].transAxes,
            color="blue",
            fontsize=fsize,
        )
        ax[i, 0].set_ylabel("Displacement (m)", fontsize=fsize)
        # ax[i, 0].set_xlim(P_arr - 150.0, S_arr + 300.0)
        ax[0, 1].legend(
            ncol=2,
            prop={"size": fsize},
            loc="center left",
            bbox_to_anchor=(0.6, 0.94),
            bbox_transform=fig.transFigure,
        )
        ax[0, 0].legend(
            ncol=2,
            prop={"size": fsize},
            loc="center left",
            bbox_to_anchor=(0.1, 0.94),
            bbox_transform=fig.transFigure,
        )
        if P_zoom is not None:
            ax[i, 0].set_ylim(-P_zoom, P_zoom)
            ax[i, 0].set_xlim(-20.0, 50.0)
        if S_zoom is not None:
            ax[i, 1].set_ylim(-S_zoom, S_zoom)
            ax[i, 1].set_xlim(-20.0, 70.0)

        ax[-1, 0].set_xlabel("Time (s)", fontsize=fsize)
        ax[-1, 1].set_xlabel("Time (s)", fontsize=fsize)
    fig.text(
        0.55,
        0.9,
        f"{event.name}\n BP:{left_BP}-{right_BP}",
        ha="right",
        va="bottom",
        size="medium",
        color="black",
        fontsize=30,
    )
    plt.savefig(pjoin(save_folder, f"{event.name}_PS.svg"))

    """ Plot P and S separate for each component """
    fig, ax = plt.subplots(nrows=3, ncols=4, sharex="col", sharey="col", figsize=(30, 15))
    Hodo_P = obspy.Stream()
    Hodo_S = obspy.Stream()
    for i, comp in enumerate(components):
        tr = st[i].copy()
        if left_BP is not None and right_BP is not None:
            _PreProcess.filter_tr(tr, fmin=left_BP, fmax=right_BP, zerophase=False)
        o_time = event.origin_time - tr.stats.starttime

        """ Plot the origin time, P and S arrivals """

        ls = "dashed"
        if P_shift != 0.0:
            ax[i, 0].axvline(x=-P_shift, c="black", ls=ls, label="P-arrival (original)")
            ls = "dotted"
        ax[i, 0].axvline(x=0, c="black", ls=ls, label="P-arrival")
        ls = "dashed"
        if S_shift != 0.0:
            ax[i, 2].axvline(x=-S_shift, c="black", ls=ls, label="S-arrival (original)")
            ls = "dotted"
        ax[i, 2].axvline(x=0, c="black", ls=ls, label="S-arrival")

        """ Slice data """
        tr_P = tr.slice(
            starttime=event.origin_time + P_arr - P_pre,
            endtime=event.origin_time + P_arr + P_post,
        )
        tr_S = tr.slice(
            starttime=event.origin_time + S_arr - S_pre,
            endtime=event.origin_time + S_arr + S_post,
        )

        """ Plot sliced data """
        # ax[i, 0].plot(
        #     tr_P.times() - P_pre, tr_P.data, label=f"selected P data", c="red",
        # )
        # ax[i, 2].plot(
        #     tr_S.times() - P_pre, tr_S.data, label=f"selected S data", c="blue",
        # )
        ax[i, 0].scatter(tr_P.times() - P_pre, tr_P.data, c=tr_P.times(), cmap="seismic")
        ax[i, 2].scatter(tr_S.times() - S_pre, tr_S.data, c=tr_S.times(), cmap="seismic")

        """ Set plotting parameters """
        fsize = 15
        ax[i, 0].text(
            s=f"P{comp}",
            x=0.99,
            y=0.90,
            ha="right",
            transform=ax[i, 0].transAxes,
            color="blue",
            fontsize=fsize,
        )
        ax[i, 2].text(
            s=f"S{comp}",
            x=0.99,
            y=0.90,
            ha="right",
            transform=ax[i, 2].transAxes,
            color="blue",
            fontsize=fsize,
        )
        ax[i, 0].set_ylabel("Displacement (m)", fontsize=fsize)
        # ax[i, 0].set_xlim(P_arr - 150.0, S_arr + 300.0)
        ax[0, 2].legend(
            ncol=2,
            prop={"size": fsize},
            loc="center left",
            bbox_to_anchor=(0.6, 0.94),
            bbox_transform=fig.transFigure,
        )
        ax[0, 0].legend(
            ncol=2,
            prop={"size": fsize},
            loc="center left",
            bbox_to_anchor=(0.1, 0.94),
            bbox_transform=fig.transFigure,
        )
        if P_zoom is not None:
            ax[i, 0].set_ylim(-P_zoom, P_zoom)
            ax[i, 0].set_xlim(-P_pre, P_post)
        if S_zoom is not None:
            ax[i, 2].set_ylim(-S_zoom, S_zoom)
            ax[i, 2].set_xlim(-S_pre, S_post)

        ax[-1, 0].set_xlabel("Time (s)", fontsize=fsize)
        ax[-1, 2].set_xlabel("Time (s)", fontsize=fsize)

        Hodo_P += tr_P
        Hodo_S += tr_S

    """ Plot Hodograms """
    ind = np.roll(np.arange(3), 1)
    for i in range(3):
        dt = Hodo_P[i].stats.delta
        """ P-hodogram """
        ntP = len(Hodo_P[i].data)
        t = np.linspace(0, dt * ntP, ntP, endpoint=False) - P_pre
        sc = ax[i, 1].scatter(Hodo_P[i].data, Hodo_P[ind[i]].data, c=t, cmap="seismic",)
        # ax.set_xlim(st_real[0].data.min(), st_real[0].data.max())
        # ax.set_ylim(st_real[1].data.min(), st_real[1].data.max())
        ax[i, 1].set_xlabel(Hodo_P[i].stats.channel, fontsize=fsize, c="blue")
        ax[i, 1].set_ylabel(Hodo_P[ind[i]].stats.channel, fontsize=fsize, c="blue")

        """ S-hodogram """
        ntS = len(Hodo_S[i].data)
        t = np.linspace(0, dt * ntS, ntS, endpoint=False) - S_pre
        sc = ax[i, 3].scatter(Hodo_S[i].data, Hodo_S[ind[i]].data, c=t, cmap="seismic",)
        # ax.set_xlim(st_real[0].data.min(), st_real[0].data.max())
        # ax.set_ylim(st_real[1].data.min(), st_real[1].data.max())
        ax[i, 3].set_xlabel(Hodo_S[i].stats.channel, fontsize=fsize, c="blue")
        ax[i, 3].set_ylabel(Hodo_S[ind[i]].stats.channel, fontsize=fsize, c="blue")

    fig.text(
        0.55,
        0.9,
        f"{event.name}\n BP:{left_BP}-{right_BP}",
        ha="right",
        va="bottom",
        size="medium",
        color="black",
        fontsize=30,
    )
    plt.savefig(pjoin(save_folder, f"{event.name}_Hodogram_PS.svg"))

