"""
Script to determine focal mechanism of the InSight station.
Important: ssh-copy-id -i .ssh/id_rsa sysop@marshost.ethz.ch (once)
"""
__author__ = "Nienke Brinkman"

import argparse
import toml

import SS_MTI

# from SS_MTI.GetData import GetData
# from SS_MTI.PullData import PullData


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
    source = toml.load(args.input_file, _dict=dict)

    ## Step 2:
    """ Get the observed data """
    OBS = GetData()
    OBS.read_inv(inv_path=source["DATA"]["inventory_filepath"])  # Inventory file
    OBS.read_cat(cat_path=source["DATA"]["catalog_filepath"])  # Catalog file
    OBS.read_events_from_cat(
        event_params=source["EVENTS"],
        cat=OBS.cat,
        inv=OBS.inv,
        local_folder=source["DATA"]["waveform_filepath"],
        host_name=source["SERVER"]["host_name"],
        user_name=source["SERVER"]["username"],
        remote_folder=source["SERVER"]["remote_folder"],
    )

    a = 1

    ## Step 3:
    """ Start inversion """

    pass
