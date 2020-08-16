"""
Interactive example to determine focal mechanism of the InSight station.
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from os import listdir as lsdir
import instaseis
import numpy as np
import matplotlib.pyplot as plt

import SS_MTI
import EventInterface
from SS_MTI import PostProcessing as _PostProcessing

save_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Trial_2"

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
        "amplitude_correction": ["PZ"],
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
        "amplitude_correction": ["PZ"],
        "t_pre": [1, 1, 1, 1, 1],
        "t_post": [17, 30, 30, 17, 30],
        "weights": [[1, 3], [1, 3], [1, 3], [1, 3], [1, 3]],
        "start_weight_len": 7.0,
        "dt": 0.05,
        "db_path": "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE",
        "npz_file": "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz",
        "ylims": [2e-9, 4e-9, 4e-9, 2e-9, 4e-9],
    },
    "S0183a": {
        "phases": ["P", "P"],
        "components": ["Z", "R"],
        "phase_corrs": [2.1, 2.1, 2.0, -0.8, 2.0],
        "tstars": [1.2, 1.2],
        "fmin": 0.2,
        "fmax": 0.4,
        "zerophase": False,
        "amplitude_correction": ["PZ"],
        "t_pre": [1, 1],
        "t_post": [30, 30],
        "weights": [[1, 3], [1, 3]],
        "start_weight_len": 7.0,
        "dt": 0.05,
        "db_path": "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE",
        "npz_file": "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz",
        "ylims": [2e-10, 4e-10],
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

""" Specify receiver """
lat_rec = 4.502384
lon_rec = 135.623447
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

""" """
depths = np.arange(5, 90, 3)
# depths = [20]

# strikes = np.arange(0, 360, 20)
# dips = np.arange(0, 91, 15)
# rakes = np.arange(-180, 180, 15)

# strikes = [255]
# dips = [55]
# rakes = [-85]

strikes = np.arange(0, 360, 5)
dips = np.arange(0, 91, 5)
rakes = np.arange(-180, 180, 5)

""" Loop over events to invert for: """
event_nr = 0
for i, v in event_input.items():
    event = events[event_nr]
    print(event.name)
    event_nr += 1
    assert event.name == i, "Dictionary and events do not iterate correct"
    if event_nr < 2:
        continue
    """ Define forward modeler """
    forward_method = "INSTASEIS"
    db_path = v["db_path"]

    mnt_folder = "/mnt/marshost/"
    if not lsdir(mnt_folder):
        print(f"{mnt_folder} is still empty, mounting now...")
        SS_MTI.DataGetter.mnt_remote_folder(
            host_ip="marshost.ethz.ch",
            host_usr="sysop",
            remote_folder="/data/",
            mnt_folder=mnt_folder,
        )

    db = instaseis.open_db(db_path)

    # SS_MTI.DataGetter.unmnt_remote_folder(mnt_folder=mnt_folder)

    npz_file = v["npz_file"]

    """ Define misfit function """
    misfit_method = "L2"

    weights = v["weights"]
    start_weight_len = v["start_weight_len"]
    dt = v["dt"]

    """ Define inversion method """
    inv_method = "GS"
    phases = v["phases"]
    components = v["components"]
    amplitude_correction = v["amplitude_correction"]
    t_pre = v["t_pre"]
    t_post = v["t_post"]
    phase_corrs = v["phase_corrs"]
    tstars = v["tstars"]
    fmin = v["fmin"]
    fmax = v["fmax"]
    zerophase = v["zerophase"]
    output_folder = save_folder
    ylims = v["ylims"]

    """ Extra phases to plot:"""
    extra_phases = None  # ["PP", "SS", "pP", "sP", "PPP", "SSS"]

    if forward_method == "INSTASEIS":
        fwd = SS_MTI.Forward.Instaseis(
            instaseis_db=db,
            taup_model=npz_file,
            or_time=event.origin_time,
            dt=dt,
            start_cut=100.0,
            end_cut=800.0,
        )
    elif forward_method == "REFLECTIVITY":
        fwd = SS_MTI.Forward.reflectivity()
    else:
        raise ValueError(
            "forward_method can be either INSTASEIS or REFLECTIVITY in [FORWARD] of .toml file"
        )

    if misfit_method == "L2":
        misfit = SS_MTI.Misfit.L2(weights=weights, start_weight_len=start_weight_len, dt=dt)
    elif misfit_method == "CC":
        misfit = SS_MTI.Misfit.CC(shift_samples=128)
    elif misfit_method == "POL":
        misfit = SS_MTI.Misfit.POL()
    else:
        raise ValueError("misfit can be L2, CC or POL in [MISFIT] of .toml file")

    """ Start inversion """
    # if inv_method == "GS":
    SS_MTI.Inversion.Grid_Search_run(
        fwd=fwd,
        misfit=misfit,
        event=event,
        rec=rec,
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
        list_to_correct_M0=amplitude_correction,
        output_folder=output_folder,
        plot=True,
        plot_extra_phases=extra_phases,
        color_plot="blue",
        Ylims=ylims,
    )
    # elif inv_method == "Direct":
    # """ Direct inversion """
    SS_MTI.Inversion.Direct(
        fwd=fwd,
        misfit=misfit,
        event=event,
        rec=rec,
        phases=phases,
        components=components,
        phase_corrs=phase_corrs,
        t_pre=t_pre,
        t_post=t_post,
        depths=depths,
        tstars=tstars,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        output_folder=output_folder,
        plot=True,
        plot_extra_phases=extra_phases,
        color_plot="red",
        Ylims=ylims,
    )
    # else:
    #     raise ValueError("inv_method is not recognized, specify: GS or Direct")

    """ Post-processing """

    """ (waveform plotting post inversion from generated files)"""
    # _PostProcessing.post_waveform_plotting(
    #     h5_file_folder=output_folder,
    #     method="GS",
    #     misfit_name=misfit.name,
    #     misfit_weight_len=misfit.start_weight_len,
    #     fwd=fwd,
    #     event=event,
    #     rec=rec,
    #     phases=phases,
    #     components=components,
    #     t_pre=t_pre,
    #     t_post=t_post,
    #     depths=[41],
    #     phase_corrs=phase_corrs,
    #     fmin=fmin,
    #     fmax=fmax,
    #     zerophase=zerophase,
    #     tstars=tstars,
    #     plot_extra_phases=extra_phases,
    #     Ylims=ylims,
    # )

    """ (misfit vs depth analysis)"""
    DOF = sum([int((x + y) / v["dt"]) for x, y in zip(v["t_pre"], v["t_post"])])
    Moho_d = 30
    fig = _PostProcessing.plot_misfit_vs_depth(
        save_paths=[output_folder],
        event_name=event.name,
        DOF=DOF,
        depths=depths,
        misfit_name=misfit.name,
        veloc_model=fwd.veloc_name,
        true_depth=None,
        Moho=Moho_d,
        fmin=fmin,
        fmax=fmax,
        amount_of_phases=len(v["phases"]),
    )
    plt.tight_layout()
    plt.savefig(
        pjoin(
            save_folder,
            f"Misfit_vs_Depth_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.svg",
        ),
        dpi=600,
    )
    plt.close()

    """ (best MT vs depth phase arrivals) """
    # depths = depths[::2]  # np.array([23, 26, 29])  #
    # t_pre = [5, 5]
    # t_post = [40, 40]
    # phases = [phases[0], phases[1]]
    # components = [components[0], components[1]]
    # phase_corrs = [phase_corrs[0], phase_corrs[1]]
    # tstars = [tstars[0], tstars[1]]
    # fig = _PostProcessing.plot_phases_vs_depth(
    #     h5_file_folder=output_folder,
    #     method="GS",
    #     misfit_name=misfit.name,
    #     fwd=fwd,
    #     event=event,
    #     rec=rec,
    #     phases=phases,
    #     components=components,
    #     t_pre=t_pre,
    #     t_post=t_post,
    #     depths=depths,
    #     phase_corrs=phase_corrs,
    #     fmin=fmin,
    #     fmax=fmax,
    #     zerophase=zerophase,
    #     tstars=tstars,
    #     color_plot="blue",
    # )
    # # plt.tight_layout()
    # plt.savefig(
    #     pjoin(
    #         save_folder,
    #         f"PhaseTracking_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.svg",
    #     ),
    #     dpi=600,
    # )
    # plt.close()

    """ Uncertainty estimates:"""
    # fig = _PostProcessing.Source_Uncertainty(
    #     h5_file_folder=output_folder,
    #     event_name=event.name,
    #     method="GS",
    #     misfit_name=misfit.name,
    #     fwd=fwd,
    #     phases=phases,
    #     components=components,
    #     depths=depths,
    #     DOF=DOF,
    #     fmin=fmin,
    #     fmax=fmax,
    # )
    # plt.tight_layout()
    # plt.savefig(
    #     pjoin(
    #         save_folder,
    #         f"Uncertainties_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.svg",
    #     ),
    #     dpi=600,
    # )
    # plt.close()

