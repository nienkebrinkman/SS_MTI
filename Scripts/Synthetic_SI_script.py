"""
Interactive example to determine focal mechanism of the InSight station.
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from obspy.taup import TauPyModel
import obspy
import instaseis

import SS_MTI
from EventInterface import EventObj


## Step 1: Define parameters
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 10.99032013
lon_src = 170
depth = 45.0
name = "Test_Event"

phases = ["S", "P"]

mnt_folder = "/mnt/marshost/"

SS_MTI.DataGetter.mnt_remote_folder(
    host_ip="marshost.ethz.ch", host_usr="sysop", remote_folder="/data/", mnt_folder=mnt_folder,
)


npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz"
# npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
model = TauPyModel(npz_file)

db_path = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
# db_path = "http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/"
db = instaseis.open_db(db_path)

SS_MTI.DataGetter.unmnt_remote_folder(mnt_folder=mnt_folder)


lat_rec = 4.502384
lon_rec = 135.623447

strike = 60
dip = 90
rake = 0
focal_mech = [strike, dip, rake]
M0 = 5.62e13

dt = 0.05

components = "ZRT"
kind = "displacement"
noise = True

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
    M0=M0,
    dt=dt,
    components=components,
    kind=kind,
    noise=noise,
)
event = event.event
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

## Step 3:
""" Define forward modeler """
fwd = SS_MTI.Forward.Instaseis(
    instaseis_db=db,
    taup_model=npz_file,
    or_time=event.origin_time,
    dt=dt,
    start_cut=100.0,
    end_cut=800.0,
)

## Step 4:
""" Define misfit """
misfit_method = "L2"

weights = [[1, 3], [1, 3]]
start_weight_len = 7.0


if misfit_method == "L2":
    misfit = SS_MTI.Misfit.L2(weights=weights, start_weight_len=start_weight_len, dt=dt)
elif misfit_method == "CC":
    misfit = SS_MTI.Misfit.CC(shift_samples=128)
elif misfit_method == "POL":
    misfit = SS_MTI.Misfit.POL()
else:
    raise ValueError("misfit can be L2, CC or POL in [MISFIT] of .toml file")

## Step 5:
""" Start inversion """
components = ["T", "Z"]
# components = ["Z"]
amplitude_correction = ["PZ", "ST"]
t_pre = [1, 1]
t_post = [20, 20]
depths = [depth]
strikes = [180, strike]
dips = [dip, 10, 20, 30]
rakes = [rake]
phase_corrs = None
tstars = None
fmin = 1.0 / 8.0
fmax = 1.0 / 5.0
zerophase = False
output_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/"

""" Grid-Search inversion """
# SS_MTI.Inversion.Grid_Search_run(
#     fwd=fwd,
#     misfit=misfit,
#     event=event,
#     rec=rec,
#     phases=phases,
#     components=components,
#     t_pre=t_pre,
#     t_post=t_post,
#     depths=depths,
#     strikes=strikes,
#     dips=dips,
#     rakes=rakes,
#     phase_corrs=phase_corrs,
#     tstars=tstars,
#     fmin=fmin,
#     fmax=fmax,
#     zerophase=zerophase,
#     list_to_correct_M0=amplitude_correction,
#     output_folder=output_folder,
#     plot=True,
# )

# """ Direct inversion """
SS_MTI.Inversion.Direct(
    fwd=fwd,
    misfit=misfit,
    event=event,
    rec=rec,
    phases=phases,
    components=components,
    t_pre=t_pre,
    t_post=t_post,
    depths=depths,
    tstars=tstars,
    fmin=fmin,
    fmax=fmax,
    zerophase=zerophase,
    output_folder=output_folder,
    plot=True,
)

