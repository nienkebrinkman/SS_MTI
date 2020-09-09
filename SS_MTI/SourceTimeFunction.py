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


def stf_tstar(tstar, dt, npts):
    # Create source time function with
    # M(f) = 1 / (1+(f/fc)**2)
    f = fftpack.fftfreq(npts, dt)
    stf_amp = np.exp(-np.pi * np.abs(f) * tstar)

    # Set phase to obtain a minimum phase signal
    # Phi is the hilbert transform of
    # the log of the spectral amplitude
    stf_phase = fftpack.hilbert(np.log(stf_amp))

    stf_td = fftpack.ifft(stf_amp * np.exp(1j * stf_phase))
    return stf_td.real, f, stf_amp


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