""" Script to compute synthetic time & frequency ASCII files """

from os.path import join as pjoin
from os import listdir as lsdir
import numpy as np
import obspy
import instaseis
from matplotlib import mlab as mlab
import matplotlib.pyplot as plt

from SS_MTI.Read_H5 import Read_GS_h5, Read_Direct_Inversion
from SS_MTI import Forward, DataGetter

"""  Parameters """
save_folder = "/home/nienke/Documents/Research/Data/MTI/Synthetic_spectra/"
folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Trial_4/"

event_name = "S0183a"
phases = ["P"]
tstar = [1.2]
depth = 29
fmin = 0.2
fmax = 0.4
misfit_name = "L2" 
amount_of_phases = 5

dt = 0.05

component = "Z"

veloc_name =  "TAYAK_BKE"
db_path = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
npz_file = f"/home/nienke/Documents/Research/Data/npz_files/{veloc_name}.npz"

""" Open file within preferred depth range """
GS_file_name = pjoin(
                folder,
                f"GS_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_name}.hdf5",
            )

depth_GS, sdr, M0_GS, misfit_L2_GS = Read_GS_h5(
    Filename=GS_file_name, amount_of_phases=amount_of_phases
)
lowest_ind = np.sum(misfit_L2_GS, axis=1).argsort()[0]
depth = depth_GS[lowest_ind]
MT = sdr[lowest_ind,:]
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
event_input = {
    event_name: []}
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
    winlen = min(winlen_sec * Fs,
                 (tr.stats.endtime -
                  tr.stats.starttime) * Fs / 2.)
    NFFT = obspy.signal.util.next_pow_2(winlen)
    pad_to = np.max((NFFT * 2, 1024))
    p, f = mlab.psd(tr.data,
                    Fs=Fs, NFFT=NFFT, detrend='linear',
                    pad_to=pad_to, noverlap=NFFT // 2)
    return f, p

""" Generate Green's functions per depth """
fig, ax = plt.subplots(nrows=len(phases), ncols=2, figsize=(8, 5 * len(phases)))
for i, phase in enumerate(phases):
    syn_GF = fwd.get_greens_functions(
        comp=component,
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
        filter=True,
        fmin=fmin,
        fmax=fmax,
        zerophase=False,
    )
    syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)

    tr = fwd.generate_synthetic_data(
                st_GF=syn_GF, focal_mech=MT, M0=M0, slice=True,tt = syn_tt,t_pre = 1, t_post=30,
            )
    time = np.vstack((tr.times()-1,tr.data))
    with open(pjoin(save_folder,f'{event.name}_{veloc_name}_{phase}_time.txt'), 'wb') as file:
        np.save(file, time, allow_pickle=False)

    win_len_sec = 20.
    f, p = calc_PSD(tr,winlen_sec=win_len_sec)
    freq = np.vstack((f,p))

    with open(pjoin(save_folder,f'{event.name}_{veloc_name}_{phase}_spectra.txt'), 'wb') as file:
        np.save(file, freq, allow_pickle=False)

    if len(phases) == 1:
        ax[0].plot(tr.times()- 1., tr.data)
        ax[0].set_xlabel("Time (s)")
        ax[0].set_ylabel("Displacement (m)")
        ax[0].axis("tight")
        ax[0].set_title(phase)

        ax[1].plot(f,p)
        ax[1].set_xlabel("Frequency (Hz)")
        ax[1].set_ylabel("Power Spectral Density")
        ax[1].axis("tight")
        ax[1].set_title(phase)

    else:
        ax[i,0].plot(tr.times()- 1., tr.data)
        ax[i,0].set_xlabel("Time (s)")
        ax[i,0].set_ylabel("Displacement (m)")
        ax[i,0].axis("tight")
        ax[i,0].set_title(phase)

        ax[i,1].plot(f,p)
        ax[i,1].set_xlabel("Frequency (Hz)")
        ax[i,1].set_ylabel("Power Spectral Density")
        ax[i,1].axis("tight")
        ax[i,1].set_title(phase)

plt.savefig(pjoin(save_folder,f"{event.name}_{veloc_name}.pdf"))

