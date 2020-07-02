"""
Interactive example to determine focal mechanism of the InSight station.
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist

import SS_MTI
import EventInterface



or_time= obspy.UTCDateTime("2020-3-10T12:00:00"),
lat= 10.99032013,
lon= 170,
depth= 45.0,
name= "Test_Event",

event = EventInterface.create_event(or_time=or_time,lat=lat,lon=lon,depth=depth,name=name)

## Variables that look if the event file is already saved from previous runs:
inv = SS_MTI.DataGetter.read_inv(inv_path = "/path/to/inv")  # Inventory file
cat = SS_MTI.DataGetter.read_cat(cat_path=f_in["DATA"]["catalog_filepath"])  # Catalog file
events = SS_MTI.DataGetter.read_events_from_cat(
    event_params=f_in["EVENTS"],
    cat=cat,
    inv=inv,
    local_folder=f_in["DATA"]["waveform_filepath"],
    host_name=f_in["SERVER"]["host_name"],
    user_name=f_in["SERVER"]["username"],
    remote_folder=f_in["SERVER"]["remote_folder"],
    save_file_name=cat_save_name,
)

## Step 3:
""" Define forward modeler """

forward_method = f_in["FORWARD"]["METHOD"]
forward_dict = f_in["FORWARD"][forward_method]

if forward_method == "INSTASEIS":
    fwd = SS_MTI.Forward.Instaseis(
        instaseis_db=forward_dict["VELOC"],
        taup_model=forward_dict["VELOC_taup"],
        rec_lat=f_in["PARAMETERS"]["RECEIVER"],
        rec_lon=rec_lon,
    )
elif forward_method == "REFLECTIVITY":
    fwd = SS_MTI.Forward.reflectivity()
else:
    raise ValueError(
        "forward_method can be either INSTASEIS or REFLECTIVITY in [FORWARD] of .toml file"
    )

fwd = SS_MTI.Forward()




""" Start inversion """


invs = SS_MTI.Inversion(
    forward_method=forward_method,
    forward_dict=forward_dict,
    rec_lat=f_in["PARAMETERS"]["RECEIVER"]["la_r"],
    rec_lon=f_in["PARAMETERS"]["RECEIVER"]["lon_r"],
)
inv_methods = f_in["INVERSION"]["METHOD"]
for inv_method in inv_methods:
    print("Start {} inversion".format(inv_method))
    for event in OBS.events:
        if inv_method == 'GS'
            invs.Grid_Search(event= event,depths=[10],strikes=[10], dips=[10], rakes=[10])
            
        
        pass
