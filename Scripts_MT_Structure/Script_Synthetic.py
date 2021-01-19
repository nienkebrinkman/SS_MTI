"""
Interactive example to invert focal mechanism and structure simultaneously
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from os import listdir as lsdir
from obspy.taup import TauPyModel
import obspy
import instaseis
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.cross_correlation import xcorr_max, correlate

import SS_MTI
from EventInterface import EventObj
from SS_MTI import PostProcessing as _PostProcessing


## Step 1: Define parameters
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 0
lon_src = 0
depth = 40.0
name = "Test_Event"

phases = ["P", "S", "S", "P", "S"]

mnt_folder = "/mnt/marshost/"

if not lsdir(mnt_folder):
    print(f"{mnt_folder} is still empty, mounting now...")
    SS_MTI.DataGetter.mnt_remote_folder(
        host_ip="marshost.ethz.ch",
        host_usr="sysop",
        remote_folder="/data/",
        mnt_folder=mnt_folder,
    )


npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
# npz_file = "/home/nienke/Data_2020/npz_files/TAYAK_BKE.npz"
# npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
model = TauPyModel(npz_file)

db_path = "/mnt/marshost/instaseis/databases/blindtestmodels_1s/TAYAK_1s"
# db_path = "/opt/databases/TAYAK_15s_BKE"
# db_path = "http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/"
db = instaseis.open_db(db_path)

# SS_MTI.DataGetter.unmnt_remote_folder(mnt_folder=mnt_folder)


lat_rec = -30
lon_rec = 0

strike = 0
dip = 45  # 45
rake = 0  # -90
M0 = 5.62e13
# [m_rr, m_tt, m_pp, m_rt, m_rp, m_tp] = SS_MTI.GreensFunctions.convert_SDR(strike, dip, rake, M0)
[m_rr, m_tt, m_pp, m_rt, m_rp, m_tp] = [0.0, 1.0e17, 0.0, 0.0, 0.0, 0.0]
focal_mech = [m_rr, m_tt, m_pp, m_rt, m_rp, m_tp]

dt = 0.025

components = "ZNE"
kind = "displacement"
noise = False

Parallel = False

## Step 2:
"""Create observed data and waveforms """
event = EventObj(
    or_time=or_time,
    lat_src=lat_src,
    lon_src=lon_src,
    lat_rec=lat_rec,
    lon_rec=lon_rec,
    depth=depth,
    name=name,
)
event.add_picks(taup_model=model, depth=depth, phases=phases)
event.add_waveforms(
    instaseis_db_path=db_path,
    focal_mech=focal_mech,
    M0=None,
    dt=dt,
    components=components,
    kind=kind,
    noise=noise,
)
event = event.event
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

## Get the P and S arrivals
P_arr = event.picks["P"] - event.origin_time
S_arr = event.picks["S"] - event.origin_time

## Plot the synthetic event
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(28, 12))
# colours = ["r", "b", "k"]
# for i, comp in enumerate(["Z", "R", "T"]):
#     ax.axvline(x=P_arr, c="red", ls="dashed", label="P-arrival")
#     ax.axvline(x=S_arr, c="blue", ls="dashed", label="S-arrival")
#     st_select = event.waveforms_VBB.select(channel="BH" + comp)
#     print(st_select[0].stats.starttime)
#     for j, tr in enumerate(st_select):
#         ax.plot(
#             tr.times(), tr.data + j * 1e-9, colours[i], label=f"Comp: {comp}",
#         )
#         # ax.plot([P_arr, P_arr], [-0.5e-9 + j * 1e-9, 0.5e-9 + j * 1e-9])
# plt.legend()
# plt.show()


# ###
# north_coor = np.array([0, 0, 1])
# src_coor = np.array(
#     [
#         np.cos(np.deg2rad(lat_src)) * np.cos(np.deg2rad(lon_src)),
#         np.cos(np.deg2rad(lat_src)) * np.sin(np.deg2rad(lon_src)),
#         np.sin(np.deg2rad(lat_src)),
#     ]
# )
# # src_coor = np.array(
# #     [
# #         np.sin(np.deg2rad(lat_src)) * np.cos(np.deg2rad(lon_src)),
# #         np.sin(np.deg2rad(lat_src)) * np.sin(np.deg2rad(lon_src)),
# #         np.cos(np.deg2rad(lat_src)),
# #     ]
# # )

# R = np.cross(north_coor, src_coor)
# R_norm = R / (np.sqrt(np.sum(np.square(R), axis=0)))

# rec_coor = np.array(
#     [
#         np.cos(np.deg2rad(lat_rec)) * np.cos(np.deg2rad(lon_rec)),
#         np.cos(np.deg2rad(lat_rec)) * np.sin(np.deg2rad(lon_rec)),
#         np.sin(np.deg2rad(lat_rec)),
#     ]
# )

# # R_new_coor = np.cross(rec_coor, R_norm) + 0
# # R_new_coor = np.cross(north_coor, rec_coor) + 0
# # rec_coor = np.array(
# #     [
# #         np.sin(np.deg2rad(lat_rec)) * np.cos(np.deg2rad(lon_rec)),
# #         np.sin(np.deg2rad(lat_rec)) * np.sin(np.deg2rad(lon_rec)),
# #         np.cos(np.deg2rad(lat_rec)),
# #     ]
# # )

# beta = np.deg2rad(lat_src)
# R_new_coor = (
#     np.cos(beta) * rec_coor
#     + np.sin(beta) * (np.cross(R_norm, rec_coor))
#     + (np.inner(R_norm, rec_coor)) * (1 - np.cos(beta)) * R_norm
# )

# baz = np.rad2deg(np.arctan2(R_new_coor[1], R_new_coor[0]))
# rot_baz = np.abs(lon_src - baz)
# epi = np.rad2deg(np.arctan2(np.sqrt(R_new_coor[0] ** 2 + R_new_coor[1] ** 2), R_new_coor[2]))


# polar_src = lat_src + 90
# polar_rec = lat_rec + 90
# az_src = lon_src
# az_rec = lon_rec

# x_src = np.cos(np.deg2rad(lat_src)) * np.cos(np.deg2rad(lon_src))
# y_src = np.cos(np.deg2rad(lat_src)) * np.sin(np.deg2rad(lon_src))
# z_src = np.sin(np.deg2rad(lat_src))

# x_rec = np.cos(np.deg2rad(lat_rec)) * np.cos(np.deg2rad(lon_rec))
# y_rec = np.cos(np.deg2rad(lat_rec)) * np.sin(np.deg2rad(lon_rec))
# z_rec = np.sin(np.deg2rad(lat_rec))

# x_dist = np.abs(x_src - x_rex)
# y_dist = np.abs(y_src - y_rex)
# distance = np.sqrt(x_dist ** 2 + y_dist ** 2)
# angle = np.rad2deg(np.arctan2(y_dist, x_dist))

# Determine back azimuth
lon_st = lon_rec
lon_ev = event.longitude
lat_ev = event.latitude
lat_st = lat_rec
x1_tan = np.sin((lon_st - lon_ev) * 180.0 / np.pi)
x2_tan = np.cos(lat_st * 180.0 / np.pi) * np.sin(lat_ev * 180.0 / np.pi) - np.sin(
    lat_st * 180.0 / np.pi
) * np.cos(lat_ev * 180.0 / np.pi) * np.cos((lon_st - lon_ev) * 180.0 / np.pi)
azi = np.arctan2(x1_tan, x2_tan) * 180.0 / np.pi
# azi += 180
if azi <= 0:
    bazi = azi + 180.0
elif azi > 0:
    bazi = azi - 180.0
print(bazi)

save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/"
import Create_Vmod


bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
Create_Vmod.create_dat_file(
    src_depth=depth,
    focal_mech=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    M0=None,
    epi=event.distance,
    baz=bazi,
    save_path=save_path,
    bm_file_path=bm_file_path,
)

# import subprocess

# subprocess.call("./crfl_sac_mars", shell=True, cwd=save_path)

## Open de files from the reflectivity code:
from os import listdir as lsdir
from os.path import isfile

st_refl = obspy.Stream()
st_files = [f for f in lsdir(save_path) if f.startswith("st") if isfile(pjoin(save_path, f))]

for st_file in st_files:
    st_temp = obspy.read(pjoin(save_path, st_file))

    st_temp[0].stats.channel = "xx" + st_file.split(".")[-1]
    distance = st_temp[0].stats.sac.gcarc
    st_temp[0].stats.distance = distance

    B = st_temp[0].stats.sac.b  # beginning time
    tstart = st_temp[0].stats.starttime  # absolute starttime
    et0 = tstart - B
    st_temp[0] = st_temp[0].trim(
        starttime=tstart - B, endtime=st_temp[0].stats.endtime, pad=True, fill_value=0.0
    )
    st_refl += st_temp

st_refl.select(channel="xxz")[0].data *= -1
st_refl.select(channel="xxt")[0].data *= -1
# st_refl.select(channel="xxr")[0].data *= -1

print(st_refl)

# st_refl.resample(20)

import scipy.signal as signal


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


Stream_object = event.waveforms_VBB.copy()
st_refl_cop = st_refl.copy()
# st_refl_cop = st_refl_cop.rotate("NE->RT", back_azimuth=180)


# ## Filter both
dt_ins = Stream_object[0].stats.delta
dt_ref = st_refl_cop[0].stats.delta
ff = np.array([0.03, 0.4])
Stream_object = taper_trace(Stream_object)
Stream_object = bandpass_filt(Stream_object, ff[0], ff[1], dt_ins)
st_refl_cop = taper_trace(st_refl_cop)
st_refl_cop = bandpass_filt(st_refl_cop, ff[0], ff[1], dt_ref)

# ## Scale both
# Stream_object, st_refl_cop = scale_traces(Stream_object, st_refl_cop)

fig, ax = plt.subplots(nrows=3, ncols=1, sharex="all", sharey="all", figsize=(15, 15))
for i, comp in enumerate(["Z", "R", "T"]):
    ax[i].axvline(x=P_arr, c="red", ls="dashed", label="P-arrival")
    ax[i].axvline(x=S_arr, c="blue", ls="dashed", label="S-arrival")
    ax[i].plot(
        Stream_object.traces[i].times(),
        Stream_object.traces[i].data,
        label=f"INSTASEIS {Stream_object[i].stats.channel}",
        c="k",
    )

    st_select = st_refl_cop.select(channel="xx" + comp)
    for j, tr in enumerate(st_select):
        ax[i].plot(
            tr.times(), tr.data * 0.4e4, label=f"REFLECTIVITY {comp}",
        )
        ax[i].set_xlim(P_arr - 10.0, S_arr + 100.0)
    # ax[i].set_ylim(-0.2, 0.2)
    ax[i].set_ylabel("Displacement (m)")
    ax[i].legend()
ax[-1].set_xlabel("Time (s)")
# ax[2].set_ylim(-0.2e-30,0.2e-30)
plt.show()

## Cross-correlate the reflectivity traces with the instaseis traces:
# 1. copy the stream:
st_refl_cop2 = st_refl_cop.copy()
st_ins_cop2 = Stream_object.copy()
# 2. cut out a window around the S-wave:
timing = st_refl_cop2[0].stats.starttime + S_arr
min_range = 5
max_range = 50
st_refl_cop2 = st_refl_cop2.trim(
    starttime=timing - min_range, endtime=timing + max_range, pad=True, fill_value=0.0
)
timing = event.origin_time + S_arr
st_ins_cop2 = st_ins_cop2.trim(
    starttime=timing - min_range, endtime=timing + max_range, pad=True, fill_value=0.0
)
fig, ax = plt.subplots(nrows=3, ncols=1, sharex="all", sharey="all", figsize=(15, 15))
for i, (comp_ins, comp_refl) in enumerate(zip(["Z", "N", "E"], ["Z", "R", "T"])):
    tr_ins = st_ins_cop2.select(channel="BH" + comp_ins).traces[0]
    tr_refl = st_refl_cop2.select(channel="xx" + comp_refl).traces[0]
    tr_refl.data *= 0.4e4

    # tr_refl.data *= 0.0
    # tr_refl.data[10:20] = 1.0

    # tr_ins.data *= 0.0
    # tr_ins.data[30:40] = 1.0

    # 3. Correlate
    corrarray = correlate(tr_refl, tr_ins, domain="time", shift=128)
    shift_CC, misfit_CC = xcorr_max(corrarray, abs_max=False)
    print(shift_CC)
    # 4. Shift
    # misfit_CC[iphase] = corrarray[(len(corrarray) - 1) // 2 + int(shifts[phases[iphase]])]
    # shift_CC[iphase] = shifts[phases[iphase]]

    ax[i].axvline(x=0, c="blue", ls="dashed", label="S-arrival")
    ax[i].plot(
        tr_ins.times() - min_range, tr_ins.data, label=f"INSTASEIS {tr_ins.stats.channel}", c="k",
    )

    ax[i].plot(
        tr_refl.times() - min_range - shift_CC * dt_ref,
        tr_refl.data,
        label=f"REFLECTIVITY {tr_refl.stats.channel}",
    )
    ax[i].plot(
        tr_refl.times() - min_range,
        tr_refl.data,
        label=f"No-shift REFLECTIVITY {tr_refl.stats.channel}",
    )
    ax[i].set_ylabel("Displacement (m)")
    ax[i].legend()
ax[-1].set_xlabel("Time (s)")
# ax[2].set_ylim(-0.2e-30,0.2e-30)
plt.show()


## Step 3:
# """ Define forward modeler """
# fwd = SS_MTI.Forward.reflectivity(
#     or_time=event.origin_time, dt=dt, start_cut=100.0, end_cut=800.0,
# )

# """ Create a synthetic seismogram """
# syn_GF = fwd.get_greens_functions(
#     comp=["Z"],
#     depth=depth,
#     distance=event.distance,
#     lat_src=event.latitude,
#     lon_src=event.longitude,
#     rec=rec,
#     tstar=None,
#     LQT=False,
#     inc=None,
#     baz=baz,
#     M0=M0,
#     filter=False,
#     fmin=None,
#     fmax=None,
#     zerophase=zerophase,
# )

# tr_syn_full = fwd.generate_synthetic_data(
#     st_GF=syn_GFs[i], focal_mech=focal_mech, M0=M0, slice=False,
# )
