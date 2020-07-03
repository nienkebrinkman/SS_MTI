"""
Script to determine focal mechanism of the InSight station.
Important: ssh-copy-id -i .ssh/id_rsa sysop@marshost.ethz.ch (once)
"""
__author__ = "Nienke Brinkman"

import argparse
import toml
from os.path import join as pjoin
from os.path import exists as exist

import SS_MTI

def define_arguments():
    helptext = "Determine focal mechanisms of Marsquake"
    parser = argparse.ArgumentParser(description=helptext)

    helptext = "Input toml file"
    parser.add_argument("input_file", help=helptext)
    return parser.parse_args()


if __name__ == "__main__":
    args = define_arguments()

    ## Step 1:
    """ Read input file """
    f_in = toml.load(args.input_file, _dict=dict)

    # If you want to save your catalog file for only the events that you want to use in your inversion:
    cat_save_name = args.input_file.split("/")[-1].strip(".toml")  # None if you dont want to save

    ## Step 2:
    """ Get the observed data """

    ## Variables that look if the event file is already saved from previous runs:
    inv = SS_MTI.DataGetter.read_inv(inv_path=f_in["DATA"]["inventory_filepath"])  # Inventory file
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

    event = event[0]

    ## Step 3:
    """ Define forward modeler """

    forward_method = f_in["FORWARD"]["METHOD"]
    forward_dict = f_in["FORWARD"][forward_method]

    if forward_method == "INSTASEIS":
        fwd = SS_MTI.Forward.Instaseis(
            instaseis_db=forward_dict["VELOC"],
            taup_model=forward_dict["VELOC_taup"],
            rec_lat=f_in["PARAMETERS"]["RECEIVER"]["la_r"],
            rec_lon=f_in["PARAMETERS"]["RECEIVER"]["lon_r"],
            or_time=event.origin_time,
            dt = event.delta
            start_cut=f_in["PARAMETERS"]["start_cut"],
            end_cut=f_in["PARAMETERS"]["end_cut"],
        )
    elif forward_method == "REFLECTIVITY":
        fwd = SS_MTI.Forward.reflectivity()
    else:
        raise ValueError(
            "forward_method can be either INSTASEIS or REFLECTIVITY in [FORWARD] of .toml file"
        )
   
    ## Step 3:
    """ Define misfit """
    misfit_method = f_in["MISFIT"]["METHOD"]
    misfit_dict = f_in["MISFIT"][misfit_method]

    if method == "L2":
        misfit = SS_MTI.Misfit.L2()
    elif method == "CC":
        misfit = SS_MTI.Misfit.CC()
    elif method == "POL":
        misfit = SS_MTI.Misfit.POL()
    else:
        raise ValueError("misfit can be L2, CC or POL in [MISFIT] of .toml file")


    ## step 4:
    """ Start inversion """

    inv_methods = f_in["INVERSION"]["METHOD"]
    for inv_method in inv_methods:
        print("Start {} inversion".format(inv_method))
        for event in OBS.events:
            if inv_method == 'GS'
                invs.Grid_Search(event= event,depths=[10],strikes=[10], dips=[10], rakes=[10])
                
         
            pass
