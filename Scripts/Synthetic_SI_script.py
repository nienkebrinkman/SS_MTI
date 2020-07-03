"""
Interactive example to determine focal mechanism of the InSight station.
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from obspy.taup import TauPyModel
import obspy

import SS_MTI
from EventInterface import EventObj

## Step 1: Define parameters
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 10.99032013
lon_src = 170
depth = 45.0
name = "Test_Event"

phases = ["P", "S"]

npz_file = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"
model = TauPyModel(npz_file)

db_path = "http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/"


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
noise = False

## Step 2: Create observed data and waveforms
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

## Step 3:
""" Define forward modeler """
fwd = SS_MTI.Forward.Instaseis(
    instaseis_db=db_path,
    taup_model=npz_file,
    rec_lat=lat_rec,
    rec_lon=lon_rec,
    or_time=event.origin_time,
    dt=dt,
    start_cut=100.0,
    end_cut=800.0,
)

## Step 4:
""" Define misfit """
misfit_method = "L2"

weights = [[1, 3], [1, 4]]
start_weight_len = 3.0


if misfit_method == "L2":
    misfit = SS_MTI.Misfit.L2(weights=weights, start_weight_len=start_weight_len, dt=dt)
elif misfit_method == "CC":
    misfit = SS_MTI.Misfit.CC()
elif misfit_method == "POL":
    misfit = SS_MTI.Misfit.POL()
else:
    raise ValueError("misfit can be L2, CC or POL in [MISFIT] of .toml file")

## Step 5:
""" Start inversion """
components = ["Z", "Z"]
t_pre = [1, 1]
t_post = [10, 10]
depths = [10, 11, 12]
strikes = [40, 50]
dips = [20, 30]
rakes = [100, -100]
phase_corrs = None
tstars = None
fmin = 1.0 / 8.0
fmax = 1.0 / 5.0
zerophase = False

SS_MTI.Inversion.Grid_Search_run(
    fwd=fwd,
    misfit=misfit,
    event=event,
    phases=phases,
    components=components,
    t_pre=t_pre,
    t_post=t_post,
    depths=depths,
    strikes=strikes,
    dips=dips,
    rakes=rakes,
    phase_corrs=phase_corrs,
    tstars=tstars,
    fmin=fmin,
    fmax=fmax,
    zerophase=zerophase,
)


# invs = SS_MTI.Inversion(
#     forward_method=forward_method,
#     forward_dict=forward_dict,
#     rec_lat=f_in["PARAMETERS"]["RECEIVER"]["la_r"],
#     rec_lon=f_in["PARAMETERS"]["RECEIVER"]["lon_r"],
# )
# inv_methods = f_in["INVERSION"]["METHOD"]
# for inv_method in inv_methods:
#     print("Start {} inversion".format(inv_method))
#     for event in OBS.events:
#         if inv_method == "GS":
#             invs.Grid_Search(event=event, depths=[10], strikes=[10], dips=[10], rakes=[10])

#             pass
