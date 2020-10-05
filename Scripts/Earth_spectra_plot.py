#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2019
:license:
    None
"""
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import obspy
import scipy.fftpack
from matplotlib import mlab as mlab
import numpy as np

from mqs_reports.catalog import Catalog


def calc_PSD(tr, winlen_sec):
    # tr.taper(0.05)
    Fs = tr.stats.sampling_rate
    winlen = min(winlen_sec * Fs, (tr.stats.endtime - tr.stats.starttime) * Fs / 2.0)
    NFFT = obspy.signal.util.next_pow_2(winlen)
    pad_to = np.max((NFFT * 2, 1024))
    p, f = mlab.psd(tr.data, Fs=Fs, NFFT=NFFT, detrend="linear", pad_to=pad_to, noverlap=NFFT // 2)
    return f, p


st_earth = obspy.read("/home/nienke/Documents/Research/Data/Earth/Package_1576834180297.mseed")
inv_earth = obspy.read_inventory(
    "/home/nienke/Documents/Research/Data/Earth/Package_1576834180297_0.xml"
)
st_earth = st_earth.select(channel="HHZ")
st_earth.remove_response(inv_earth, output="DISP")

coords = [37.20, 22.80, 69.0]
from obspy.geodetics import locations2degrees

tr_obs_filt = st_earth.select(channel="HHZ")[1]
tr_obs_raw = tr_obs_filt.copy()
fmin = 0.1  # 1.0 / 30.0
fmax = 2.0
tr_obs_filt.filter("highpass", freq=fmin)
tr_obs_filt.filter("highpass", freq=fmin)
tr_obs_filt.filter("lowpass", freq=fmax)

arrs_P = [183.5, 243, 197.0, 170.0]
P = arrs_P[2]
phases = ["P"]
t_pres = [60]
t_posts = [400]
Normalize = False
win_len_sec = [10.0]

# slice the traces:
tr_obs_filt = tr_obs_filt.slice(
    starttime=tr_obs_filt.stats.starttime + P,
    endtime=tr_obs_filt.stats.starttime + P + t_posts[0],
)
tr_obs_raw = tr_obs_raw.slice(
    starttime=tr_obs_raw.stats.starttime + P, endtime=tr_obs_filt.stats.starttime + P + t_posts[0],
)

fig, ax = plt.subplots(nrows=len(phases), ncols=2, figsize=(18, 6 * len(phases)))


f_obs_filt, p_obs_filt = calc_PSD(tr_obs_filt, winlen_sec=win_len_sec[0])
f_obs_raw, p_obs_raw = calc_PSD(tr_obs_raw, winlen_sec=win_len_sec[0])

if Normalize:
    # tr_syn_filt.normalize()
    # tr_syn_raw.normalize()
    tr_obs_raw.normalize()
    tr_obs_filt.normalize()

if len(phases) == 1:
    max_val = np.abs(tr_obs_raw.max())
    max_val = max_val + 0.2 * max_val
    if Normalize:
        max_val = 2.0

    # ax[0].plot(
    #     tr_syn_raw.times() - t_pres[i],
    #     tr_syn_raw.data,
    #     color="red",
    #     label="raw synthetic",
    #     alpha=0.7,
    # )
    # ax[0].plot(
    #     tr_syn_filt.times() - t_pres[i],
    #     tr_syn_filt.data + max_val,
    #     color="red",
    #     ls="--",
    #     label="filtered synthetic",
    # )
    ax[0].plot(
        tr_obs_raw.times() - t_pres[0],
        tr_obs_raw.data,
        color="black",
        label="raw observed",
        alpha=0.7,
    )
    ax[0].plot(
        tr_obs_filt.times() - t_pres[0],
        tr_obs_filt.data + max_val,
        color="black",
        ls="--",
        label="filtered observed",
    )
    ax[0].set_xlabel("Time (s)", fontsize=25)

    if Normalize:
        ax[0].set_ylabel("Displacement", fontsize=25)
    else:
        ax[0].set_ylabel("Displacement (m)", fontsize=25)
    ax[0].axis("tight")
    # ax[0].set_title(f"{phase}{components[i]}", color="blue")
    ax[0].get_yaxis().get_offset_text().set_visible(False)
    ax_max = max(ax[0].get_yticks())
    exponent_axis = np.floor(np.log10(ax_max)).astype(int)
    # ax[0].annotate(
    #     r"$\times$10$^{%i}$" % (exponent_axis),
    #     xy=(0.01, 0.9),
    #     xycoords="axes fraction",
    #     fontsize=20,
    # )
    ax[0].text(
        s=r"$\times$10$^{%i}$" % (exponent_axis),
        x=0.02,
        y=0.93,
        ha="left",
        transform=ax[0].transAxes,
        color="black",
        fontsize=20,
    )
    ymax = ax[0].get_ylim()[0]
    ax[0].axvline(0, color="blue")
    ax[0].text(
        0 - 1, ymax * 0.85, "P", verticalalignment="center", color="blue", fontsize=25,
    )

    # ax[0].axvline(S, color="blue")
    # ax[0].text(
    #     S - 1, ymax * 0.85, "S", verticalalignment="center", color="blue", fontsize=25,
    # )
    if Normalize:
        ax[0].axes.get_yaxis().set_ticks([])

    # ax[1].semilogy(f_syn_raw, p_syn_raw, color="red", alpha=0.7, label="raw synthetic")
    # ax[1].semilogy(f_syn_filt, p_syn_filt, color="red", ls="--", label="filtered synthetic")
    ax[1].semilogx(
        f_obs_raw, 10 * np.log10(p_obs_raw), color="black", alpha=0.7, label="raw observed"
    )
    ax[1].semilogx(
        f_obs_filt, 10 * np.log10(p_obs_filt), color="black", ls="--", label="filtered observed"
    )
    ax[1].axvline(fmin, color="black")
    ax[1].axvline(fmax, color="black")
    ax[1].axvspan(fmin, fmax, facecolor="orange", alpha=0.3)
    ax[1].set_xlabel("Frequency (Hz)", fontsize=25)
    ax[1].set_ylabel("displacement PSD [dB]", fontsize=25)
    ax[1].axis("tight")
    # ax[1].set_title(f"{phase}{components[i]} t* used: {tstar[i]}", color="blue")
    ax[1].set_xlim(5e-2, 5e0)
    # ax[1].set_ylim(1e-24, 1e-16)
    ax[1].set_ylim(-250, -100)

    ax[0].text(
        s="Greece",
        x=0.98,
        y=0.9,
        ha="right",
        transform=ax[0].transAxes,
        color="blue",
        fontsize=25,
    )

    ax[0].tick_params(axis="both", which="major", labelsize=18)
    ax[1].tick_params(axis="both", which="major", labelsize=18)
    ax[0].tick_params(axis="both", which="minor", labelsize=10)
    ax[1].tick_params(axis="both", which="minor", labelsize=10)

    ax[1].legend(fontsize=20)

plt.savefig(
    "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_2/Spectra_paperfig/Earth.svg",
    dpi=600,
)

