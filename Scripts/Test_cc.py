from obspy.signal.cross_correlation import xcorr_max, correlate
import numpy as np
import matplotlib.pyplot as plt

a = np.array([1, 1, 1, 1, 1, 3e14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
b = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

corrarray = correlate(a, b, domain="time", shift=20)
shift_CC, misfit_CC = xcorr_max(corrarray, abs_max=False)

c = 1

