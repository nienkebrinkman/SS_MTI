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


save_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Sigma_est_weight_noise_new"

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
        "phase_corrs": [0.2, 10.1, 10.7, 0.2, 10.7],
        "tstars": [0.8, 1.1, 1.1, 0.8, 1.1],
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
        "phase_corrs": [-0.5, 2.5, 2.5, -0.5, 2.5, -0.5],
        "tstars": [1.1, 1.2, 1.2, 1.1, 1.2, 1.1],
        "fmin": 0.1,
        "fmax": 0.7,
        "zerophase": False,
        "amplitude_correction": ["PZ"],
        "t_pre": [1, 1, 1, 1, 1, 1],
        "t_post": [17, 30, 30, 17, 30, 17],
        "weights": [[1, 3], [1, 3], [1, 3], [1, 3], [1, 3], [1, 3]],
        "start_weight_len": 7.0,
        "dt": 0.05,
        "db_path": "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE",
        "npz_file": "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz",
        "ylims": [2e-9, 4e-9, 4e-9, 2e-9, 4e-9, 2e-9],
    },
    "S0183a": {
        "phases": ["P", "P"],
        "components": ["Z", "R"],
        "phase_corrs": [2.1, 2.1],
        "tstars": [1.5, 1.5],
        "fmin": 1.0 / 5.0,
        "fmax": 1.0 / 2.5,
        "zerophase": False,
        "amplitude_correction": ["PZ", "PR"],
        "t_pre": [1, 1],
        "t_post": [30, 30],
        "weights": [[1, 3], [1, 3]],
        "start_weight_len": 7.0,
        "dt": 0.05,
        "db_path": "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE",
        "npz_file": "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz",
        "ylims": [2e-10, 2e-10],
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
lat_rec = 4.5  # 02384
lon_rec = 135.623447
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

""" """
depths = np.arange(5, 90, 3)
# depths = np.arange(29, 50, 3)
# depths = [8]

strikes = np.arange(0, 360, 20)
dips = np.arange(0, 91, 15)
rakes = np.arange(-180, 180, 15)

# strikes = [15.0116557194]  # [132.395557582]
# dips = [59.551091053]  # [51.9591191063]
# rakes = [-45.6275510954]  # [-139.94976385]

# strikes = np.arange(0, 360, 5)
# dips = np.arange(0, 91, 5)
# rakes = np.arange(-180, 180, 5)

""" Define different velocity models"""
db_name_1 = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
npz_file_name_1 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz"

db_name_2 = "/mnt/marshost/instaseis2/databases/TAYAK_shallow"
npz_file_name_2 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"

db_names = [db_name_1]  # , db_name_3, db_name_4, db_name_5]
npz_file_names = [npz_file_name_1]

""" Loop over events to invert for: """
event_nr = 0
for i, v in event_input.items():
    event = events[event_nr]
    print(event.name)
    event_nr += 1
    assert event.name == i, "Dictionary and events do not iterate correct"
    if event.name == "S0173a":
        pass
    else:
        continue

    if event.name == "S0183a":
        event.distance = 44.5
        event.baz = 73
        event.az = 253
        event.latitude = 15.09
        event.longitude = 179.59

    """ Define forward modeler """
    forward_method = "INSTASEIS"
    # db_path = v["db_path"]
    # npz_file = v["npz_file"]
    db_nr = 0
    for db_path, npz_file in zip(db_names, npz_file_names):
        db_nr += 0
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
        extra_phases = None  # [
        #     "PP",
        #     "sP",
        #     "pP",
        # ]

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
            misfit = SS_MTI.Misfit.Pol(
                components=components, start_weight_len=start_weight_len, weights=weights, dt=dt
            )
        else:
            raise ValueError("misfit can be L2, CC or POL in [MISFIT] of .toml file")

        """ Start inversion """

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
        #     depths=depths,
        #     phase_corrs=phase_corrs,
        #     fmin=fmin,
        #     fmax=fmax,
        #     zerophase=zerophase,
        #     tstars=tstars,
        #     plot_extra_phases=extra_phases,
        #     Ylims=ylims,
        # )

        # _PostProcessing.post_waveform_plotting(
        #     h5_file_folder=output_folder,
        #     method="Direct",
        #     misfit_name=misfit.name,
        #     misfit_weight_len=misfit.start_weight_len,
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
        #     plot_extra_phases=extra_phases,
        #     Ylims=ylims,
        # )

        # """ (misfit vs depth analysis)"""
        DOF = sum([int((x + y) / v["dt"]) for x, y in zip(v["t_pre"], v["t_post"])])
        Moho_d = 24
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
        # depths_phases = depths[1::2]  # np.array([23, 26, 29])  #
        # t_pre = [5, 5]
        # t_post = [40, 65]
        # phases = [phases[0], phases[1]]
        # components = [components[0], components[1]]
        # phase_corrs = [phase_corrs[0], phase_corrs[1]]
        # tstars = [tstars[0], tstars[1]]
        # # tstars = [tstar_P, tstar_S]
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
        #     depths=depths_phases,
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
        # fig, fig_sdr = _PostProcessing.Source_Uncertainty(
        #     h5_file_folder=output_folder,
        #     event_name=event.name,
        #     method="GS",
        #     misfit_name=misfit.name,
        #     fwd=fwd,
        #     phases=phases,
        #     components=components,
        #     depths=np.arange(26, 47, 3),
        #     DOF=DOF,
        #     fmin=fmin,
        #     fmax=fmax,
        # )
        # fig.tight_layout()
        # fig.savefig(
        #     pjoin(
        #         save_folder,
        #         f"Uncertainties_FULL_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.svg",
        #     ),
        #     dpi=600,
        # )
        # plt.close(fig)
        # fig_sdr.tight_layout()
        # fig_sdr.savefig(
        #     pjoin(
        #         save_folder,
        #         f"Uncertainties_SDR_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.svg",
        #     ),
        #     dpi=600,
        # )
        # plt.close(fig_sdr)

