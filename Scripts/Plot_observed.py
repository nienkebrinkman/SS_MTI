"""
Plot observed data from Marsquakes and an Earthquake
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
import instaseis
import numpy as np

import SS_MTI
import EventInterface

save_folder = "/home/nienke/Documents/Research/Data/MTI/Observed"

path = "/home/nienke/Documents/Research/Data/MTI/old_catalog"
# path = "/home/nienke/Documents/Research/SS_MTI/Data"
path_to_inventory = pjoin(path, "inventory.xml")
path_to_catalog = pjoin(path, "catalog.xml")


""" Read the inventory and catalog file (the once that contain info about the marsquakes) """
inv = SS_MTI.DataGetter.read_inv(inv_path=path_to_inventory)  # Inventory file
cat = SS_MTI.DataGetter.read_cat(cat_path=path_to_catalog)  # Catalog file


""" Define events to invert for and its parameters """
event_input = {
    "S0235b": {
        "phases": ["P", "S", "S", "P", "S"],
        "components": ["Z", "T", "Z", "R", "R"],
        "phase_corrs": [0.0, 10.1, 10.1, 0.0, 10.10],
        "tstars": [1.2, 1.5, 1.5, 1.2, 1.5],
        "fmin": 0.1,
        "fmax": 0.9,
        "zerophase": False,
        "amplitude_correction": ["PZ", "ST"],
        "t_pre": [1, 1, 1, 1, 1],
        "t_post": [30, 30, 30, 30, 30],
        "weights": [[1, 3], [1, 3], [1, 3], [1, 3], [1, 3]],
        "start_weight_len": 7.0,
        "dt": 0.05,
        "db_path": "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE",
        "npz_file": "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz",
        "ylims": [4e-9, 4e-9, 4e-9, 4e-9, 4e-9],
    },
    "S0173a": {
        "phases": ["P", "S", "S", "P", "S"],
        "components": ["Z", "T", "Z", "R", "R"],
        "phase_corrs": [-0.8, 2.0, 2.0, -0.8, 2.0],
        "tstars": [1.2, 1.6, 1.6, 1.2, 1.6],
        "fmin": 0.1,
        "fmax": 0.7,
        "zerophase": False,
        "amplitude_correction": ["PZ", "ST"],
        "t_pre": [1, 1, 1, 1, 1],
        "t_post": [17, 30, 30, 17, 30],
        "weights": [[1, 3], [1, 3], [1, 3], [1, 3], [1, 3]],
        "start_weight_len": 7.0,
        "dt": 0.05,
        "db_path": "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE",
        "npz_file": "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz",
        "ylims": [2e-9, 4e-9, 4e-9, 2e-9, 4e-9],
    },
}

""" Get the data into a list of obspy.Event objects """
events = SS_MTI.DataGetter.read_events_from_cat(
    event_params=event_input,
    cat=cat,
    inv=inv,
    local_folder="/mnt/marshost/",
    host_name="marshost.ethz.ch",
    user_name="sysop",
    remote_folder="/data/",
    save_file_name=pjoin(save_folder, "event.mseed"),
)

event = events[1]

obs_tt = []
for i, phase in enumerate(phases):
    obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
st_obs, sigmas = _PreProcess.prepare_event_data(
    event=event,
    phases=phases,
    components=components,
    slice=True,
    tts=obs_tt,
    t_pre=t_pre,
    t_post=t_post_new,
    filter=filter_par,
    fmin=fmin,
    fmax=0.5,
    zerophase=zerophase,
    noise_level=False,
)

""" Specify receiver """
lat_rec = 4.502384
lon_rec = 135.623447
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

a = 1
