""" Script to compute synthetic time & frequency ASCII files """

from os.path import join as pjoin
from os import listdir as lsdir
import numpy as np
import obspy
import instaseis
from matplotlib import mlab as mlab
import matplotlib.pyplot as plt
from obspy import UTCDateTime as utct

from SS_MTI.Read_H5 import Read_GS_h5, Read_Direct_Inversion
from SS_MTI import Forward, DataGetter
from SS_MTI import PreProcess as _PreProcess

"""  Parameters """
save_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_1/Synthetic_spectra/"
folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_1/5phases_weightchange/"

# event_name = "S0235b"
# phases = ["P", "S", "S", "P", "S"]
# phase_corrs = [0.2, 10.1, 10.1, 0.2, 10.7]
# components = ["Z", "T", "Z", "R", "R"]
# tstar = [0.4, 0.4, 0.4, 0.4, 0.4]
# t_pres = [1, 1, 1, 1, 1]
# t_posts = [30, 30, 30, 30, 30]
# depth = 29
# fmin = 0.1
# fmax = 0.9
# misfit_name = "L2"
# amount_of_phases = 5

event_name = "S0173a"
phases = ["P", "S", "S", "P", "S"]
phase_corrs = [-0.5, 2.5, 1.5, -0.5, 2.5]
components = ["Z", "T", "Z", "R", "R"]
tstar = [0.4, 0.4, 0.4, 0.4, 0.4]
t_pres = [1, 1, 1, 1, 1]
t_posts = [30, 30, 30, 30, 30]
depth = 29
fmin = 0.1
fmax = 0.7
misfit_name = "L2"
amount_of_phases = 5

dt = 0.05

veloc_name = "TAYAK_BKE"
db_path = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
npz_file = f"/home/nienke/Documents/Research/Data/npz_files/{veloc_name}.npz"

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
st_obs, sigmas_noise = _PreProcess.prepare_event_data(
    event=event,
    phases=phases,
    components=components,
    slice=True,
    tts=obs_tt,
    t_pre=t_pres,
    t_post=t_posts,
    filter=True,
    fmin=fmin,
    fmax=fmax,
    zerophase=False,
    noise_level=True,
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
    Fs = tr.stats.sampling_rate
    winlen = min(winlen_sec * Fs, (tr.stats.endtime - tr.stats.starttime) * Fs / 2.0)
    NFFT = obspy.signal.util.next_pow_2(winlen)
    pad_to = np.max((NFFT * 2, 1024))
    p, f = mlab.psd(tr.data, Fs=Fs, NFFT=NFFT, detrend="linear", pad_to=pad_to, noverlap=NFFT // 2)
    return f, p


""" Generate Green's functions per depth """
fig, ax = plt.subplots(nrows=len(phases), ncols=2, figsize=(8, 5 * len(phases)))
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
        M0=1.0,
        filter=False,
        fmin=fmin,
        fmax=fmax,
        zerophase=False,
    )
    syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)

    tr = fwd.generate_synthetic_data(
        st_GF=syn_GF,
        focal_mech=MT,
        M0=M0,
        slice=True,
        tt=syn_tt,
        t_pre=t_pres[i],
        t_post=t_posts[i],
    )
    time = np.vstack((tr.times() - 1, tr.data))
    with open(pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}_time.txt"), "wb") as file:
        np.save(file, time, allow_pickle=False)

    win_len_sec = 20.0
    f, p = calc_PSD(tr, winlen_sec=win_len_sec)
    freq = np.vstack((f, p))

    with open(pjoin(save_folder, f"{event.name}_{veloc_name}_{phase}_spectra.txt"), "wb") as file:
        np.save(file, freq, allow_pickle=False)

    f_obs, p_obs = calc_PSD(st_obs[i], winlen_sec=win_len_sec)

    if len(phases) == 1:
        ax[0].plot(tr.times() - t_pres[i], tr.data)
        ax[0].plot(st_obs[i].times() - t_pres[i], st_obs[i].data)
        ax[0].set_xlabel("Time (s)")
        ax[0].set_ylabel("Displacement (m)")
        ax[0].axis("tight")
        ax[0].set_title(f"{phase}{components[i]}")

        ax[1].plot(f, p)
        ax[1].plot(f_obs, p_obs)
        ax[1].set_xlabel("Frequency (Hz)")
        ax[1].set_ylabel("Power Spectral Density")
        ax[1].axis("tight")
        ax[1].set_title(f"{phase}{components[i]}")
        ax[1].set_xlim(0, 1)
        ax[1].yscale("log")

    else:
        ax[i, 0].plot(tr.times() - t_pres[i], tr.data)
        ax[i, 0].plot(st_obs[i].times() - t_pres[i], st_obs[i].data)
        ax[i, 0].set_xlabel("Time (s)")
        ax[i, 0].set_ylabel("Displacement (m)")
        ax[i, 0].axis("tight")
        ax[i, 0].set_title(f"{phase}{components[i]}")

        ax[i, 1].plot(f, p)
        ax[i, 1].plot(f_obs, p_obs)
        ax[i, 1].set_xlabel("Frequency (Hz)")
        ax[i, 1].set_ylabel("Power Spectral Density")
        ax[i, 1].axis("tight")
        ax[i, 1].set_title(f"{phase}{components[i]}")
        ax[i, 1].set_xlim(0, 1)
        ax[i, 1].yscale("log")


plt.savefig(pjoin(save_folder, f"{event.name}_{veloc_name}.pdf"))

