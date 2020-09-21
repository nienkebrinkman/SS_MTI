""" Script to compute synthetic time & frequency ASCII files """

from os.path import join as pjoin
from os import listdir as lsdir
import numpy as np
import obspy
import instaseis
from matplotlib import mlab as mlab
import matplotlib.pyplot as plt
from obspy import UTCDateTime as utct
import scipy.fftpack

from SS_MTI.Read_H5 import Read_GS_h5, Read_Direct_Inversion
from SS_MTI import Forward, DataGetter
from SS_MTI import PreProcess as _PreProcess

"""  Parameters """
save_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_1/Synthetic_spectra/"
folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_1/5phases_weightchange/"

# event_name = "S0235b"
# phases = ["P", "S", "S", "P", "S"]
# # phase_corrs = [0.2, 10.5, 11.1, 0.2, 11.1]
# phase_corrs = [0.2, 10.1, 10.7, 0.2, 10.7]
# components = ["Z", "T", "Z", "R", "R"]
# tstar = [0.4, 0.2, 0.2, 0.4, 0.2]
# # tstar = [None, None, None, None, None]
# t_pres = [1, 1, 1, 1, 1]
# t_posts = [30, 30, 30, 30, 30]
# depth = 62
# fmin = 0.1
# fmax = 0.9
# misfit_name = "L2"
# amount_of_phases = 5

event_name = "S0173a"
phases = ["P", "S", "S", "P", "S"]
phase_corrs = [-0.3, 2.9, 2.0, -0.3, 2.9]
components = ["Z", "T", "Z", "R", "R"]
tstar = [
    0.3,
    0.2,
    0.2,
    0.3,
    0.2,
]
t_pres = [1, 1, 1, 1, 1]
t_posts = [17, 30, 30, 17, 30]
depth = 50
fmin = 0.1
fmax = 0.5
misfit_name = "L2"
amount_of_phases = 5

filter_par = True
zerophase = False
win_len_sec = [10.0, 10.0, 10.0, 10.0, 10.0]
dt = 0.05

veloc_name = "TAYAK_BKE"
db_path = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
# db_path = "/mnt/marshost/instaseis2/databases/TAYAK_1s_30km"
# db_path = "/mnt/marshost/instaseis2/databases/EH45TcoldCrust1b"
# db_path = "/mnt/marshost/instaseis2/databases/TAYAK_shallow"
npz_file = f"/home/nienke/Documents/Research/Data/npz_files/{veloc_name}.npz"
# npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_30km.npz"
# npz_file = "/home/nienke/Documents/Research/Data/npz_files/EH45TcoldCrust1b.npz"

""" Open file within preferred depth range """
GS_file_name = pjoin(
    folder, f"GS_{event_name}_{depth}_{fmin}_{0.7}_{misfit_name}_{veloc_name}.hdf5",
)

depth_GS, sdr, M0_GS, misfit_L2_GS = Read_GS_h5(
    Filename=GS_file_name, amount_of_phases=amount_of_phases
)
lowest_ind = np.sum(misfit_L2_GS, axis=1).argsort()[0]
depth = depth_GS[lowest_ind]
MT = sdr[lowest_ind, :]
M0 = M0_GS[lowest_ind]

""" Mount to Mars-server """
mnt_folder = "/mnt/marshost/"
if not lsdir(mnt_folder):
    print(f"{mnt_folder} is still empty, mounting now...")
    DataGetter.mnt_remote_folder(
        host_ip="marshost.ethz.ch",
        host_usr="sysop",
        remote_folder="/data/",
        mnt_folder=mnt_folder,
    )


""" Get event data """
path = "/home/nienke/Documents/Research/Data/MTI/old_catalog"
# path = "/home/nienke/Documents/Research/SS_MTI/Data"
path_to_inventory = pjoin(path, "inventory.xml")
path_to_catalog = pjoin(path, "catalog.xml")
inv = DataGetter.read_inv(inv_path=path_to_inventory)  # Inventory file
cat = DataGetter.read_cat(cat_path=path_to_catalog)  # Catalog file
event_input = {event_name: []}
event = DataGetter.read_events_from_cat(
    event_params=event_input,
    cat=cat,
    inv=inv,
    local_folder="/mnt/marshost/",
    host_name="marshost.ethz.ch",
    user_name="sysop",
    remote_folder="/data/",
    save_file_name=pjoin(folder, "event.mseed"),
)[0]

obs_tt = []
for i, phase in enumerate(phases):
    obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
st_obs_full, sigmas_noise = _PreProcess.prepare_event_data(
    event=event,
    phases=phases,
    components=components,
    slice=False,
    tts=obs_tt,
    filter=False,
    noise_level=True,
)

st_obs_filt = obspy.Stream()
st_obs_raw = obspy.Stream()
for i, tr in enumerate(st_obs_full):
    tr_copy = tr.copy()
    _PreProcess.filter_tr(tr_copy, fmin=fmin, fmax=fmax, zerophase=zerophase)
    tr_window = tr_copy.slice(
        starttime=event.origin_time + obs_tt[i] - t_pres[i],
        endtime=event.origin_time + obs_tt[i] + t_posts[i],
    )
    st_obs_filt += tr_window

    st_obs_raw += tr.slice(
        starttime=event.origin_time + obs_tt[i] - t_pres[i],
        endtime=event.origin_time + obs_tt[i] + t_posts[i],
    )


""" Specify receiver """
lat_rec = 4.5  # 02384
lon_rec = 135.623447
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

"""Get the forward object """
db = instaseis.open_db(db_path)

fwd = Forward.Instaseis(
    instaseis_db=db,
    taup_model=npz_file,
    or_time=event.origin_time,
    dt=dt,
    start_cut=100.0,
    end_cut=800.0,
)

""" Calculate spectrum"""


def calc_PSD(tr, winlen_sec):
    # tr.taper(0.05)
    Fs = tr.stats.sampling_rate
    winlen = min(winlen_sec * Fs, (tr.stats.endtime - tr.stats.starttime) * Fs / 2.0)
    NFFT = obspy.signal.util.next_pow_2(winlen)
    pad_to = np.max((NFFT * 2, 1024))
    p, f = mlab.psd(tr.data, Fs=Fs, NFFT=NFFT, detrend="linear", pad_to=pad_to, noverlap=NFFT // 2)
    return f, p


def calc_freq(tr):
    # winlen = 20.
    tr.taper(10)
    N = len(tr.data)  # number of samples
    T = tr.stats.delta
    yf = scipy.fftpack.fft(tr.data)
    p = 2.0 / N * np.abs(yf[: N // 2])
    # NFFT = obspy.signal.util.next_pow_2(winlen)
    # from obspy.signal.freqattributes import spectrum as sp
    # p = sp(tr.data,winlen,NFFT)
    f = np.linspace(0.0, 1.0 / (2.0 * T), N // 2)
    return f, p


""" Generate Green's functions per depth """
fig, ax = plt.subplots(nrows=len(phases), ncols=2, figsize=(10, 5 * len(phases)))
for i, phase in enumerate(phases):
    syn_GF = fwd.get_greens_functions(
        comp=components[i],
        depth=depth,
        distance=event.distance,
        lat_src=event.latitude,
        lon_src=event.longitude,
        rec=rec,
        tstar=tstar[i],
        LQT=False,
        inc=None,
        baz=None,
        M0=1e14,
        filter=False,
    )
    syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)

    tr_full = fwd.generate_synthetic_data(st_GF=syn_GF, focal_mech=MT, M0=M0 / 1e14, slice=False,)

    tr_copy = tr_full.copy()
    _PreProcess.filter_tr(tr_copy, fmin=fmin, fmax=fmax, zerophase=zerophase)
    tr_syn_filt = tr_copy.slice(
        starttime=event.origin_time + syn_tt - t_pres[i],
        endtime=event.origin_time + syn_tt + t_posts[i],
    )
    tr_syn_raw = tr_full.slice(
        starttime=event.origin_time + syn_tt - t_pres[i],
        endtime=event.origin_time + syn_tt + t_posts[i],
    )
    # time = np.vstack((tr.times() - 1, tr.data))
    # with open(pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}_time.txt"), "wb") as file:
    #     np.save(file, time, allow_pickle=False)

    f_syn_filt, p_syn_filt = calc_PSD(tr_syn_filt, winlen_sec=win_len_sec[i])
    f_syn_raw, p_syn_raw = calc_PSD(tr_syn_raw, winlen_sec=win_len_sec[i])

    # f_new, p_new = calc_freq(tr)

    # freq = np.vstack((f, p))

    # with open(pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}_spectra.txt"), "wb") as file:
    #     np.save(file, freq, allow_pickle=False)

    f_obs_filt, p_obs_filt = calc_PSD(st_obs_filt[i], winlen_sec=win_len_sec[i])
    f_obs_raw, p_obs_raw = calc_PSD(st_obs_raw[i], winlen_sec=win_len_sec[i])
    # f_obs_new, p_obs_new = calc_freq(st_obs[i])

    if len(phases) == 1:
        ax[0].plot(tr_syn_raw.times() - t_pres[i], tr_syn_raw.data, color="red")
        ax[0].plot(tr_syn_raw.times() - t_pres[i], tr_syn_raw.data, color="red", ls="--")
        ax[0].plot(st_obs_raw[i].times() - t_pres[i], st_obs_raw[i].data, color="black")
        ax[0].plot(st_obs_filt[i].times() - t_pres[i], st_obs_filt[i].data, color="black", ls="--")
        ax[0].set_xlabel("Time (s)")
        ax[0].set_ylabel("Displacement (m)")
        ax[0].axis("tight")
        ax[0].set_title(f"{phase}{components[i]}")

        ax[1].semilogy(f_syn_raw, p_syn_raw, color="red")
        ax[1].semilogy(f_syn_filt, p_syn_filt, color="red", ls="--")
        ax[1].semilogy(f_obs_raw, p_obs_raw, color="black")
        ax[1].semilogy(f_obs_filt, p_obs_filt, color="black", ls="--")
        for t in [0.1, 0.3, 0.5, 0.7, 0.9]:
            p_tstar = np.exp(-np.pi * f_syn_raw * t * 2)
            ax[1].semilogy(f_syn_raw, p_tstar, label=f"t*:{t}")
        ax[1].set_xlabel("Frequency (Hz)")
        ax[1].set_ylabel("Power Spectral Density")
        ax[1].axis("tight")
        ax[1].set_title(f"{phase}{components[i]}")
        ax[1].set_xlim(0, 1.5)
        ax[1].set_ylim(1e-24, 1e-16)

    else:
        for t in [0.1, 0.5, 1.0, 2.0]:
            p_tstar = p_syn_raw[0] * np.exp(-np.pi * f_syn_raw * t * 2)
            ax[i, 1].semilogy(f_syn_raw, p_tstar, label=f"t*:{t}", lw=3, alpha=0.2)
        max_val = max(np.abs(tr_syn_raw.max()), np.abs(st_obs_raw[i].max()))

        ax[i, 0].plot(
            tr_syn_raw.times() - t_pres[i],
            tr_syn_raw.data,
            color="red",
            label="raw synthetic",
            alpha=0.7,
        )
        ax[i, 0].plot(
            tr_syn_filt.times() - t_pres[i],
            tr_syn_filt.data + max_val,
            color="red",
            ls="--",
            label="filtered synthetic",
        )
        ax[i, 0].plot(
            st_obs_raw[i].times() - t_pres[i],
            st_obs_raw[i].data,
            color="black",
            label="raw observed",
            alpha=0.7,
        )
        ax[i, 0].plot(
            st_obs_filt[i].times() - t_pres[i],
            st_obs_filt[i].data + max_val,
            color="black",
            ls="--",
            label="filtered observed",
        )
        ax[i, 0].set_xlabel("Time (s)")
        ax[i, 0].set_ylabel("Displacement (m)")
        ax[i, 0].axis("tight")
        ax[i, 0].set_title(f"{phase}{components[i]}", color="blue")

        ax[i, 1].semilogy(f_syn_raw, p_syn_raw, color="red", alpha=0.7)
        ax[i, 1].semilogy(f_syn_filt, p_syn_filt, color="red", ls="--")
        ax[i, 1].semilogy(f_obs_raw, p_obs_raw, color="black", alpha=0.7)
        ax[i, 1].semilogy(f_obs_filt, p_obs_filt, color="black", ls="--")
        ax[i, 1].axvline(fmin, ls=":", color="darkgreen")
        ax[i, 1].axvline(fmax, ls=":", color="darkgreen")
        ax[i, 1].set_xlabel("Frequency (Hz)")
        ax[i, 1].set_ylabel("Power Spectral Density")
        ax[i, 1].axis("tight")
        ax[i, 1].set_title(f"{phase}{components[i]} t* used: {tstar[i]}", color="blue")
        ax[i, 1].set_xlim(0, 1.5)
        ax[i, 1].set_ylim(1e-24, 1e-16)

        if i == 0:
            ax[i, 0].legend()
            ax[i, 1].legend()


plt.savefig(pjoin(save_folder, f"{event.name}_{veloc_name}.pdf"))

