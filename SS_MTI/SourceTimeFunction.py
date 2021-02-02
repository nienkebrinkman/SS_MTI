#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2020
:license:
    None
"""

import numpy as np
import matplotlib.pyplot as plt
import obspy
from scipy import fftpack
import pylab


def stf_tstar(tstar, dt, npts, nfft=None):
    # Create source time function with
    # M(f) = 1 / (1+(f/fc)**2)
    f = fftpack.fftfreq(npts, dt)
    stf_amp = np.exp(-np.pi * np.abs(f) * tstar)

    # stf_amp = np.exp(-fftpack.rfftfreq(nfft) * np.pi * tstar)
    # stf_phase = fftpack.hilbert(-1j * fftpack.rfftfreq(nfft) * 2.0 * np.pi * tstar)

    # Set phase to obtain a minimum phase signal
    # Phi is the hilbert transform of
    # the log of the spectral amplitude
    stf_phase = fftpack.hilbert(-np.pi * np.abs(f) * tstar)

    stf_td = fftpack.ifft(stf_amp * np.exp(1j * stf_phase))
    # stf_td = fftpack.irfft(stf_amp * 2 * np.pi)
    # f = 1
    return stf_td.real, f, stf_amp


def convolve_stf(stf: np.array, data: np.array, nfft=None):
    # conv = np.convolve(stf, data, mode="same")
    from scipy import signal

    # recovered, remainder = signal.deconvolve(data, stf)

    from obspy.signal.util import _npts2nfft
    import math

    if nfft is None:
        nfft = _npts2nfft(npts=len(data))
    stf_conv_f = np.fft.rfft(stf, n=nfft)

    # Apply a 5 percent, at least 5 samples taper at the end.
    # The first sample is guaranteed to be zero in any case.
    tlen = max(int(math.ceil(0.05 * len(data))), 5)
    taper = np.ones_like(data)
    taper[-tlen:] = signal.hann(tlen * 2)[tlen:]
    dataf = np.fft.rfft(taper * data, n=nfft)

    f = stf_conv_f
    recovered = np.fft.irfft(dataf * f)[:nfft]
    return recovered


def Create_stf_from_file(Filepath, desired_dt):
    """ File should be same format as from scardec (i.e. only usable for known Earth events)"""

    from scipy import signal

    file = np.loadtxt(Filepath, skiprows=2)
    time = file[:, 0]
    zero_time = np.where(time == 0.0)[0][0]
    time = time[zero_time:]
    dt = time[1] - time[0]
    moment_rate = file[zero_time:, 1]

    npts = int(time[-1] / desired_dt)
    moment_rate_new = signal.resample(moment_rate, npts)
    time_new = np.linspace(time[0], time[-1], npts, endpoint=False)

    # plt.close()
    # plt.figure()
    # plt.plot(time, moment_rate,'r')
    # plt.plot(time_new, moment_rate_new, 'g')
    # plt.show()
    return moment_rate_new


if __name__ == "__main__":
    dt = 1 / 10.0
    Fs = 1 / dt
    winlen_sec = 30.0
    winlen = winlen_sec * Fs
    import obspy.signal.util as UTIL

    NFFT = UTIL.next_pow_2(winlen)

    t = np.arange(0, 80, dt)
    fig, ax = plt.subplots(1, 2)

    params = {
        "legend.fontsize": "x-large",
        "figure.figsize": (15, 15),
        "axes.labelsize": 25,
        "axes.titlesize": "x-large",
        "xtick.labelsize": 25,
        "ytick.labelsize": 25,
    }
    pylab.rcParams.update(params)

    phase_name = ["P", "S"]

    for i, tstar in enumerate([0.2, 1.0, 2.0]):
        stf, f_stf, p_stf = stf_tstar(tstar, dt=dt, npts=len(t), nfft=NFFT)
        ax[0].plot(t, stf, label="T* %.2f" % (tstar))
        # ax[0].plot(t, stf, label="T* %s: %.2f" %(phase_name[i],tstar))
        # ax[1].plot(f_stf, p_stf, 'o', label="T* %s: %.2f" %(phase_name[i],tstar))
        ax[1].plot(f_stf, p_stf, "o", label="T*: %.2f" % (tstar))
    ax[0].legend()
    ax[0].set_xlabel("Time (s)")
    ax[0].set_ylabel("Amplitude factor")
    ax[1].set_xlabel("Frequency (Hz)")

    ax[1].set_xlim(0, 5)
    ax[1].set_yscale("log")
    plt.show()
