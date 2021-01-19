#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==========================================================================
# Plot reflecitivity and instaseis seismograms
# ==========================================================================

import numpy as np
import os, sys
import instaseis
import scipy.signal as signal
import obspy
from obspy import read_inventory, read
import matplotlib.pyplot as plt

# ==========================================================================


def taper_trace(stream, max_perc=0.03, ttype="hann"):
    """
    taper traces in stream
    """
    for tr in stream:
        tr.taper(max_perc, ttype)

    return stream


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


# ==========================================================================
# INSTASEIS

db = instaseis.open_db("http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/")

# Event parameters
lat_ev = 0.0  # latitude
lon_ev = 0.0  # longitude
depth_ev = 40.1e3  # depth [m]
time_ev = obspy.UTCDateTime(2019, 1, 2, 0, 0, 0)

# Station parameters
lat_st = -30.0
lon_st = 0.0

# Moment tensor
mt_rr = 0.0  # 1.0e17
mt_tt = 1.0e17
mt_pp = 0.0  # 1.0e17
mt_rt = 0.0
mt_rp = 0.0
mt_tp = 0.0

# Retrieve the waveform
receiver = instaseis.Receiver(latitude=lat_st, longitude=lon_st, network="XB", station="ELYSE")

source = instaseis.Source(
    latitude=lat_ev,
    longitude=lon_ev,
    depth_in_m=depth_ev,
    m_rr=mt_rr,
    m_tt=mt_tt,
    m_pp=mt_pp,
    m_rt=mt_rt,
    m_rp=mt_rp,
    m_tp=mt_tp,
    origin_time=time_ev,
)

st_ins1 = db.get_seismograms(source=source, receiver=receiver)

# Determine back azimuth and rotate
x1_tan = np.sin((lon_st - lon_ev) * 180.0 / np.pi)
x2_tan = np.cos(lat_ev * 180.0 / np.pi) * np.sin(lat_st * 180.0 / np.pi) - np.sin(
    lat_ev * 180.0 / np.pi
) * np.cos(lat_st * 180.0 / np.pi) * np.cos((lon_ev - lon_st) * 180.0 / np.pi)
azi = np.arctan2(x1_tan, x2_tan) * 180.0 / np.pi
if azi <= 0:
    bazi = azi + 180.0
elif azi > 0:
    bazi = azi - 180.0

st_ins = st_ins1.rotate("NE->RT", back_azimuth=bazi)
# st_ins = st_ins1.copy()

"""
comp    = []
for tr in st_ins1:
    comp.append(str(tr.stats.channel[-1]))

# Rotate to ZRT
Z    = st_ins1[comp.index('Z')].data
N    = st_ins1[comp.index('N')].data
E    = st_ins1[comp.index('E')].data
a    = (baz_ev-180.)/360.*2.*np.pi
R    =     np.cos(a)*N  + np.sin(a)*E
T    = -1.*np.sin(a)*N  + np.cos(a)*E

st_ins  = st_ins1.copy()
st_ins[0].data = Z
st_ins[0].stats.channel = 'Z'
st_ins[1].data = R
st_ins[1].stats.channel = 'R'
st_ins[2].data = T
st_ins[2].stats.channel = 'T'
"""

dt_ins = st_ins[0].stats.delta

# ==========================================================================
# REFLECTIVITY
# Read traces and info
path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_BKE/mxx/"
st_ref = read(path + "st001.z")
st_ref += read(path + "st001.r")
st_ref += read(path + "st001.t")

st_ref[0].stats.channel = "Z"
st_ref[1].stats.channel = "R"
st_ref[2].stats.channel = "T"

ntr = len(st_ref)  # number of traces in stream
GCARC = st_ref[0].stats.sac.gcarc  # epicentral distance in degrees
B2 = st_ref[0].stats.sac.b  # beginning time
tstart2 = st_ref[0].stats.starttime  # absolute starttime
et02 = tstart2 - B2  # absolute time of event
dt_ref = st_ref[0].stats.delta

# Z component in reflecitvity is reversed
st_ref[0].data *= -1.0

# ==========================================================================
# Filter - apply same filter to both to compare them
ff = np.array([0.03, 0.4])
st_ref = taper_trace(st_ref)
st_ref = bandpass_filt(st_ref, ff[0], ff[1], dt_ref)
st_ins = taper_trace(st_ins)
st_ins = bandpass_filt(st_ins, ff[0], ff[1], dt_ins)


def scale_traces(st1, st2):
    """
    scale each trace in st2 to the maximum value of each trace in st1
    """
    for i in range(len(st1)):
        max_val = np.max(st1[i].data)
        if np.max(st2[i].data) != 0:
            st2[i].data = st2[i].data / np.max(st2[i].data) * max_val

    return st1, st2


st_ins, st_ref = scale_traces(st_ins, st_ref)

# ==========================================================================
# Plot stuff
fig = plt.figure()
fig.set_size_inches(10, 5)

axs = [fig.add_subplot(3, 1, iax + 1) for iax in range(len(st_ins))]
for itr in range(ntr):
    # Normalize both traces to make them comparable (or use scale_traces)
    # st_ins[itr].normalize()
    # st_ref[itr].normalize()
    axs[itr].plot(st_ins[itr].times(), st_ins[itr].data, color="k", linewidth=1, label="Instaseis")
    print(st_ins[itr].stats.channel)
    axs[itr].plot(
        st_ref[itr].times(reftime=et02),
        st_ref[itr].data,
        color="tab:red",
        linewidth=1,
        alpha=0.8,
        label="Reflectivity",
    )
    print(st_ref[itr].stats.channel)
    axs[itr].grid(True, linestyle="dotted")
    axs[itr].legend(loc=1, framealpha=0.8)
    axs[itr].set_ylabel(st_ins[itr].stats.channel[-1], fontsize=14)
    axs[itr].yaxis.set_label_position("right")
    axs[itr].set_xlim([200, 600])
    if itr != ntr - 1:
        plt.setp(axs[itr].get_xticklabels(), visible=False)
        axs[itr].get_shared_x_axes().join(axs[itr], axs[itr + 1])

plt.setp(axs[-1].get_xticklabels(), visible=True)
axs[-1].set_xlabel("Time [s]")

fig.subplots_adjust(hspace=0)
plt.show()

