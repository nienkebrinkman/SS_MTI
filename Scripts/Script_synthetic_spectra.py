""" Script to compute synthetic time & frequency ASCII files """
__author__ = "Nienke Brinkman"

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
save_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_2/Spectra_paperfig/"
folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_2/5phases_cluster/Test_2020/"

Normalize = False


# event_name = "S0183a"
# phases = ["P"]
# phase_corrs = [2.1, 2.1]
# components = ["Z"]
# tstar = [1.5, 1.5]
# # tstar = [None, None, None, None, None]
# t_pres = [60]
# t_posts = [400]
# depth = 29
# fmin = 0.2
# fmax = 0.4
# misfit_name = "L2"
# amount_of_phases = 2


event_name = "S0235b"
phases = ["P"]
# phase_corrs = [0.2, 10.4, 11.1, 0.2, 11.1]
phase_corrs = [0.2, 10.4, 11.1, 0.2, 11.1]
components = ["Z"]
tstar = [0.4, 0.2]
# tstar = [None, None, None, None, None]
t_pres = [60]
t_posts = [400]
depth = 41
fmin = 0.1
fmax = 0.5
misfit_name = "L2"
amount_of_phases = 5

# event_name = "S0173a"
# phases = ["P"]
# phase_corrs = [-0.3, 2.9]
# components = ["Z"]
# tstar = [
#     0.1,
#     0.2,
#     0.2,
#     0.3,
#     0.2,
# ]
# t_pres = [60]
# t_posts = [400]
# depth = 38
# fmin = 0.1
# fmax = 0.4
# misfit_name = "L2"
# amount_of_phases = 5

# event_name = "S0235b"
# phases = ["P", "S", "S", "P", "S"]
# # phase_corrs = [0.2, 10.4, 11.1, 0.2, 11.1]
# phase_corrs = [0.2, 10.4, 11.1, 0.2, 11.1]
# components = ["Z", "T", "Z", "R", "R"]
# tstar = [0.4, 0.2, 0.2, 0.4, 0.2]
# # tstar = [None, None, None, None, None]
# t_pres = [1, 1, 1, 1, 1]
# t_posts = [30, 30, 30, 30, 30]
# depth = 59
# fmin = 0.1
# fmax = 0.5
# misfit_name = "L2"
# amount_of_phases = 5

# event_name = "S0173a"
# phases = ["P", "S", "S", "P", "S"]
# phase_corrs = [-0.3, 2.9, 2.0, -0.3, 2.9]
# components = ["Z", "T", "Z", "R", "R"]
# tstar = [
#     0.3,
#     0.2,
#     0.2,
#     0.3,
#     0.2,
# ]
# t_pres = [1, 1, 1, 1, 1]
# t_posts = [17, 30, 30, 17, 30]
# depth = 32
# fmin = 0.1
# fmax = 0.4
# misfit_name = "L2"
# amount_of_phases = 5

filter_par = True
zerophase = False
win_len_sec = [20.0, 20.0, 20.0, 20.0, 20.0]
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
    folder, f"GS_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_name}.hdf5",
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


S = utct(event.picks["S"]) + phase_corrs[1] - utct(event.picks["P"])

st_obs_filt = obspy.Stream()
st_obs_raw = obspy.Stream()
st_obs_pre_noise = obspy.Stream()
for i, tr in enumerate(st_obs_full):
    tr_copy = tr.copy()
    _PreProcess.filter_tr(tr_copy, fmin=fmin, fmax=fmax, zerophase=zerophase)
    tr_window = tr_copy.slice(
        starttime=event.origin_time + obs_tt[i] - t_pres[i],
        endtime=event.origin_time + obs_tt[i] + t_posts[i],
    )
    st_obs_filt += tr_window

    tr_pre_noise = tr.copy()
    st_obs_pre_noise += tr_pre_noise.slice(
        starttime=event.origin_time - 100.0, endtime=event.origin_time - 10.0
    )

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
fig, ax = plt.subplots(nrows=len(phases), ncols=2, figsize=(26, 6 * len(phases)))
for i, phase in enumerate(phases):
    syn_GF = fwd.get_greens_functions(
        comp=components[i],
        depth=depth,
        distance=event.distance,
        lat_src=event.latitude,
        lon_src=event.longitude,
        rec=rec,
        tstar=tstar[1],
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
    td_syn_filt = np.vstack((tr_syn_filt.times() - 1, tr_syn_filt.data))
    with open(
        pjoin(
            save_folder, f"{event.name}_{veloc_name}_{phase}{components[i]}_syn_Filtered_td.txt"
        ),
        "wb",
    ) as file:
        np.save(file, td_syn_filt, allow_pickle=False)

    td_syn_raw = np.vstack((tr_syn_raw.times() - 1, tr_syn_raw.data))
    with open(
        pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}{components[i]}_syn_Raw_td.txt"),
        "wb",
    ) as file:
        np.save(file, td_syn_raw, allow_pickle=False)

    if event_name == "S0183a":
        f_syn_filt, p_syn_filt = calc_PSD(tr_syn_filt, winlen_sec=win_len_sec[i])
        f_syn_raw, p_syn_raw = calc_PSD(tr_syn_raw, winlen_sec=win_len_sec[i])
    else:
        t_pre_window = 10
        t_post_window = 50
        syn_tt_window = fwd.get_phase_tt(phase="S", depth=depth, distance=event.distance)
        # if event_name == "S0173a":
        #     t_post_window = 17
        tr_syn_filt_freq = tr_syn_filt.copy()
        tr_syn_filt_freq = tr_syn_filt_freq.slice(
            starttime=event.origin_time + syn_tt_window - t_pre_window,
            endtime=event.origin_time + syn_tt_window + t_post_window,
        )
        tr_syn_raw_freq = tr_syn_raw.copy()
        tr_syn_raw_freq = tr_syn_raw_freq.slice(
            starttime=event.origin_time + syn_tt_window - t_pre_window,
            endtime=event.origin_time + syn_tt_window + t_post_window,
        )
        f_syn_filt, p_syn_filt = calc_PSD(tr_syn_filt_freq, winlen_sec=win_len_sec[i])
        f_syn_raw, p_syn_raw = calc_PSD(tr_syn_raw_freq, winlen_sec=win_len_sec[i])

    fd_syn_filt = np.vstack((f_syn_filt, p_syn_filt))
    with open(
        pjoin(
            save_folder, f"{event.name}_{veloc_name}_{phase}{components[i]}_syn_Filtered_fd.txt"
        ),
        "wb",
    ) as file:
        np.save(file, fd_syn_filt, allow_pickle=False)

    fd_syn_raw = np.vstack((f_syn_raw, p_syn_raw))
    with open(
        pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}{components[i]}_syn_Raw_fd.txt"),
        "wb",
    ) as file:
        np.save(file, fd_syn_raw, allow_pickle=False)

    # freq = np.vstack((f, p))
    # with open(pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}_spectra.txt"), "wb") as file:
    #     np.save(file, freq, allow_pickle=False)
    if event_name == "S0183a":
        f_obs_filt, p_obs_filt = calc_PSD(st_obs_filt[i], winlen_sec=win_len_sec[i])
        f_obs_raw, p_obs_raw = calc_PSD(st_obs_raw[i], winlen_sec=win_len_sec[i])
        f_obs_pre_noise, p_obs_pre_noise = calc_PSD(st_obs_pre_noise[i], winlen_sec=win_len_sec[i])

    else:
        t_pre_window = 10
        t_post_window = 50
        obs_tt_window = utct(event.picks["S"]) - event.origin_time + phase_corrs[1]
        # if event_name == "S0173a":
        #     t_post_window = 17
        st_obs_filt_freq = st_obs_filt.copy()
        st_obs_filt_freq[i] = st_obs_filt_freq[i].slice(
            starttime=event.origin_time + obs_tt_window - t_pre_window,
            endtime=event.origin_time + obs_tt_window + t_post_window,
        )
        st_obs_raw_freq = st_obs_raw.copy()
        st_obs_raw_freq[i] = st_obs_raw_freq[i].slice(
            starttime=event.origin_time + obs_tt_window - t_pre_window,
            endtime=event.origin_time + obs_tt_window + t_post_window,
        )

        f_obs_filt, p_obs_filt = calc_PSD(st_obs_filt_freq[i], winlen_sec=win_len_sec[i])
        f_obs_raw, p_obs_raw = calc_PSD(st_obs_raw_freq[i], winlen_sec=win_len_sec[i])
        f_obs_pre_noise, p_obs_pre_noise = calc_PSD(st_obs_pre_noise[i], winlen_sec=win_len_sec[i])

    # f_obs_new, p_obs_new = calc_freq(st_obs[i])

    if Normalize:
        tr_syn_filt.normalize()
        tr_syn_raw.normalize()
        st_obs_raw[i].normalize()
        st_obs_filt[i].normalize()

    if len(phases) == 1:
        cols = ["dodgerblue", "gold", "limegreen", "pink"]
        for t_i, t in enumerate([0.1, 0.5, 1.0, 2.0]):

            if event_name == "S0173a":
                p_tstar = p_syn_raw[2] * np.exp(-np.pi * f_syn_raw * t * 2)
            else:
                p_tstar = p_syn_raw[0] * np.exp(-np.pi * f_syn_raw * t * 2)
            ax[1].semilogx(
                f_syn_raw,
                10 * np.log10(p_tstar),
                c=cols[t_i],
                ls="--",
                label=f"t*:{t}",
                lw=2,
                alpha=0.6,
            )

        max_val = max(np.abs(tr_syn_raw.max()), np.abs(st_obs_raw[i].max()))
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
            st_obs_raw[i].times() - t_pres[i],
            st_obs_raw[i].data,
            color="black",
            label="raw observed",
            lw=2,
        )
        ax[0].plot(
            st_obs_filt[i].times() - t_pres[i],
            st_obs_filt[i].data + max_val,
            color="black",
            lw=1,
            label="filtered observed",
            alpha=0.7,
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
            fontsize=25,
        )
        ymax = ax[0].get_ylim()[0]
        ax[0].axvline(0, color="blue")
        ax[0].text(
            0 - 1, ymax * 0.85, "P", verticalalignment="center", color="blue", fontsize=30,
        )

        ax[0].axvline(S, color="blue")
        ax[0].text(
            S - 1, ymax * 0.85, "S", verticalalignment="center", color="blue", fontsize=30,
        )
        if Normalize:
            ax[0].axes.get_yaxis().set_ticks([])
        print(S)
        ax[1].semilogx(
            f_syn_raw, 10 * np.log10(p_syn_raw), color="red", lw=3, label="raw synthetic"
        )
        # ax[1].semilogx(
        #     f_syn_filt, 10 * np.log10(p_syn_filt), color="red", lw=2, label="filtered synthetic"
        # )
        ax[1].semilogx(
            f_obs_raw, 10 * np.log10(p_obs_raw), color="black", lw=3, label="raw observed"
        )
        # ax[1].semilogx(
        #     f_obs_filt, 10 * np.log10(p_obs_filt), color="black", lw=2, label="filtered observed",
        # )
        ax[1].semilogx(
            f_obs_filt, 10 * np.log10(p_obs_pre_noise), lw=3, color="slateblue", label="noise",
        )
        ax[1].axvline(fmin, color="black")
        ax[1].axvline(fmax, color="black")
        ax[1].axvspan(fmin, fmax, facecolor="orange", alpha=0.3)
        ax[1].set_xlabel("Frequency (Hz)", fontsize=30)
        ax[1].set_ylabel("displacement PSD [dB]", fontsize=30)
        ax[1].axis("tight")
        # ax[1].set_title(f"{phase}{components[i]} t* used: {tstar[i]}", color="blue")
        # ax[1].set_xlim(0, 1.5)
        ax[1].set_xlim(9e-2, 4e0)
        # ax[1].set_ylim(1e-24, 1e-16)
        ax[1].set_ylim(-240, -170)

        ax[0].text(
            s=f"{event.name}",
            x=0.98,
            y=0.9,
            ha="right",
            transform=ax[0].transAxes,
            color="blue",
            fontsize=30,
        )

        ax[0].tick_params(axis="both", which="major", labelsize=26)
        ax[1].tick_params(axis="both", which="major", labelsize=26)
        ax[0].tick_params(axis="both", which="minor", labelsize=15)
        ax[1].tick_params(axis="both", which="minor", labelsize=15)

        if i == 0:
            # ax[0].legend()
            ax[1].legend(fontsize=20, ncol=2)
            # ax[1].legend(fontsize=25)

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
        ax[i, 1].semilogy(f_syn_filt, p_syn_filt, color="red", lw=2)
        ax[i, 1].semilogy(f_obs_raw, p_obs_raw, color="black", alpha=0.7)
        ax[i, 1].semilogy(f_obs_filt, p_obs_filt, color="black", lw=2)
        ax[i, 1].axvline(fmin, ls=":", color="darkgreen")
        ax[i, 1].axvline(fmax, ls=":", color="darkgreen")
        ax[i, 1].axvspan(fmin, fmax, facecolor="orange", alpha=0.3)
        ax[i, 1].set_xlabel("Frequency (Hz)")
        ax[i, 1].set_ylabel("Power Spectral Density")
        ax[i, 1].axis("tight")
        ax[i, 1].set_title(f"{phase}{components[i]} t* used: {tstar[i]}", color="blue")
        ax[i, 1].set_xlim(0, 1.5)
        ax[i, 1].set_ylim(1e-24, 1e-16)

        if i == 0:
            ax[i, 0].legend()
            ax[i, 1].legend()


plt.savefig(pjoin(save_folder, f"{event.name}_{veloc_name}.svg"), dpi=600)

