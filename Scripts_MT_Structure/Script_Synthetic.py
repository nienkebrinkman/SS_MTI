"""
Interactive example to invert focal mechanism and structure simultaneously
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from os.path import isfile
from os import listdir as lsdir
from obspy.taup import TauPyModel
import obspy
import instaseis
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.cross_correlation import xcorr_max, correlate
import scipy.signal as signal


import SS_MTI
from EventInterface import EventObj
from SS_MTI import PostProcessing as _PostProcessing


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


## Step 1: Define parameters
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 0
lon_src = 0
depth = 20.0
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


lat_rec = -20
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

## Copy event data
st_ins = event.waveforms_VBB.copy()

## Rotate event
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
# st_ins = st_ins.rotate("NE->RT", back_azimuth=event.baz)

## Filter the event
dt_ins = st_ins[0].stats.delta
ff = np.array([0.03, 0.4])
st_ins = taper_trace(st_ins)
st_ins = bandpass_filt(st_ins, ff[0], ff[1], dt_ins)

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


save_path = (
    "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_2/",
    "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/Test_4/",
)

# save_path = [
#     "/home/nienke/Documents/Research/SS_MTI/External_packages/reflectivity_Mars/SRC/test/",
#     "/home/nienke/Documents/Research/SS_MTI/External_packages/reflectivity_Mars/SRC/test/Test_1/",
#     "/home/nienke/Documents/Research/SS_MTI/External_packages/reflectivity_Mars/SRC/test/Test_2/",
#     "/home/nienke/Documents/Research/SS_MTI/External_packages/reflectivity_Mars/SRC/test/Test_3/",
#     "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_BKE/mxx/",
# ]
# import Create_Vmod


# bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
# Create_Vmod.create_dat_file(
#     src_depth=depth,
#     focal_mech=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
#     M0=None,
#     epi=event.distance,
#     baz=bazi,
#     save_path=save_path,
#     bm_file_path=bm_file_path,
# )

# import subprocess

# subprocess.call("./crfl_sac_mars", shell=True, cwd=save_path)

## Open de files from the reflectivity code:


st_refls = []
for folder in save_path:
    st_refl = obspy.Stream()
    st_files = [f for f in lsdir(folder) if f.startswith("st") if isfile(pjoin(folder, f))]

    for st_file in st_files:
        st_temp = obspy.read(pjoin(folder, st_file))

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
    st_refls.append(st_refl)

print(st_refls)

st_refl_cop = st_refl.copy()
# st_refl_cop = st_refl_cop.rotate("NE->RT", back_azimuth=180)

fig, ax = plt.subplots(nrows=3, ncols=1, sharex="all", sharey="all", figsize=(15, 15))
for j, st in enumerate(st_refls):
    st_refl = st.copy()
    dt_ref = st_refl[0].stats.delta
    ff = np.array([0.03, 0.4])
    st_refl = taper_trace(st_refl)
    st_refl = bandpass_filt(st_refl, ff[0], ff[1], dt_ref)

    # ## Scale both
    # st_ins, st_refl = scale_traces(st_ins, st_refl)

    for i, comp in enumerate(["Z", "R", "T"]):
        if j == 0:
            ax[i].axvline(x=P_arr, c="red", ls="dashed", label="P-arrival")
            ax[i].axvline(x=S_arr, c="blue", ls="dashed", label="S-arrival")
            ax[i].plot(
                st_ins[i].times(),
                st_ins[i].data,
                label=f"INSTASEIS {st_ins[i].stats.channel}",
                c="k",
            )
            ax[i].set_ylabel("Displacement (m)")

        st_select = st_refl.select(channel="xx" + comp)
        for k, tr in enumerate(st_select):
            ax[i].plot(
                tr.times(), tr.data * 0.4e4, label=f"REFLECTIVITY {st_select[k].stats.channel}",
            )
            ax[i].set_xlim(P_arr - 10.0, S_arr + 100.0)
        ax[i].legend()
        # ax[i].set_ylim(-0.2, 0.2)

    ax[-1].set_xlabel("Time (s)")
    # ax[2].set_ylim(-0.2e-30,0.2e-30)
plt.show()

# ## Cross-correlate the reflectivity traces with the instaseis traces:
# # 1. copy the stream:
# st_refl_cop2 = st_refl_cop.copy()
# st_ins_cop2 = Stream_object.copy()
# # 2. cut out a window around the S-wave:
# timing = st_refl_cop2[0].stats.starttime + S_arr
# min_range = 5
# max_range = 50
# st_refl_cop2 = st_refl_cop2.trim(
#     starttime=timing - min_range, endtime=timing + max_range, pad=True, fill_value=0.0
# )
# timing = event.origin_time + S_arr
# st_ins_cop2 = st_ins_cop2.trim(
#     starttime=timing - min_range, endtime=timing + max_range, pad=True, fill_value=0.0
# )
# fig, ax = plt.subplots(nrows=3, ncols=1, sharex="all", sharey="all", figsize=(15, 15))
# for i, (comp_ins, comp_refl) in enumerate(zip(["Z", "R", "T"], ["Z", "R", "T"])):
#     tr_ins = st_ins_cop2.select(channel="BH" + comp_ins).traces[0]
#     tr_refl = st_refl_cop2.select(channel="xx" + comp_refl).traces[0]
#     tr_refl.data *= 0.4e4

#     # tr_refl.data *= 0.0
#     # tr_refl.data[10:20] = 1.0

#     # tr_ins.data *= 0.0
#     # tr_ins.data[30:40] = 1.0

#     # 3. Correlate
#     corrarray = correlate(tr_refl, tr_ins, domain="time", shift=128)
#     shift_CC, misfit_CC = xcorr_max(corrarray, abs_max=False)
#     print(shift_CC)
#     # 4. Shift
#     # misfit_CC[iphase] = corrarray[(len(corrarray) - 1) // 2 + int(shifts[phases[iphase]])]
#     # shift_CC[iphase] = shifts[phases[iphase]]

#     ax[i].axvline(x=0, c="blue", ls="dashed", label="S-arrival")
#     ax[i].plot(
#         tr_ins.times() - min_range, tr_ins.data, label=f"INSTASEIS {tr_ins.stats.channel}", c="k",
#     )

#     ax[i].plot(
#         tr_refl.times() - min_range - shift_CC * dt_ref,
#         tr_refl.data,
#         label=f"REFLECTIVITY {tr_refl.stats.channel}",
#     )
#     ax[i].plot(
#         tr_refl.times() - min_range,
#         tr_refl.data,
#         label=f"No-shift REFLECTIVITY {tr_refl.stats.channel}",
#     )
#     ax[i].set_ylabel("Displacement (m)")
#     ax[i].legend()
# ax[-1].set_xlabel("Time (s)")
# # ax[2].set_ylim(-0.2e-30,0.2e-30)
# plt.show()


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
