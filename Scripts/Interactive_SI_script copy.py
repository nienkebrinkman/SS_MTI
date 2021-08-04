"""
Interactive example to determine focal mechanism of the InSight station.
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
from os.path import exists as exist
from os import listdir as lsdir
from os import makedirs
import instaseis
import numpy as np
import matplotlib.pyplot as plt
import toml

import SS_MTI
import EventInterface
from SS_MTI import PostProcessing as _PostProcessing

# save_folder = "/home/nienke/Documents/Research/Data/MTI/Data_2021/Event_325/5phases/"
# save_folder = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/S0235b/"
# save_folder = "/home/nienke/Documents/Research/Data/MTI/S0809a_-0.7/"
save_folder = "/home/nienke/Documents/Research/Data/MTI/S0820a/5phase_longwindow_highweight_SR/"
# save_folder = "/home/nienke/Documents/Research/Data/MTI/S0235b_Alice_compare/"
# save_folder = (
# "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_2/5phases_cluster/Test_2020/"
# )

""" Read the inventory and catalog file (the once that contain info about the marsquakes) """
# path = "/home/nienke/Documents/Research/Data/MTI/Catalog/old_catalog"
path = "/home/nienke/Documents/Research/Data/MTI/Catalog"
path_to_inventory = pjoin(path, "inventory.xml")
path_to_catalog = pjoin(path, "catalog.xml")
inv = SS_MTI.DataGetter.read_inv(inv_path=path_to_inventory)  # Inventory file
cat = SS_MTI.DataGetter.read_cat(cat_path=path_to_catalog)  # Catalog file

""" Specify input file """
# input_file = "/home/nienke/Documents/Research/SS_MTI/Input/325_5phases.toml"
# input_file = "/home/nienke/Documents/Research/SS_MTI/Input/TAYAK_BKE_tstar_update.toml"
# input_file = "/home/nienke/Documents/Research/SS_MTI/Input/S0235b.toml"
# input_file = "/home/nienke/Documents/Research/SS_MTI/Input/S0820a.toml"
# input_file = "/home/nienke/Documents/Research/SS_MTI/Input/S0809a.toml"
input_file = pjoin(save_folder, "S0820a.toml")
event_input = toml.load(input_file, _dict=dict)

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

""" Specify depth range """
depths = np.arange(5, 90, 3)
# depths = [18, 33]

strikes = np.arange(0, 360, 10)
dips = np.arange(0, 91, 10)
rakes = np.arange(-180, 180, 20)

# strikes = [270]
# dips = [50]
# rakes = [140]

# bazs = np.arange(160, 360, 20)
bazs = [140]

""" Define different velocity models"""
db_name_1 = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
npz_file_name_1 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz"

db_name_2 = "/mnt/marshost/instaseis/databases/blindtestmodels_1s/TAYAK_1s"
npz_file_name_2 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"

db_name_3 = "/mnt/marshost/instaseis2/databases/TAYAK_1s_30km"
npz_file_name_3 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_30km.npz"

db_names = [db_name_1]  # , db_name_3, db_name_4, db_name_5]
npz_file_names = [npz_file_name_1]

""" Loop over events to invert for: """
event_nr = 0
for i, v in event_input.items():
    event = events[event_nr]
    """ Event location """
    # fig = _PostProcessing.Plot_event_location(
    #     la_s=event.latitude,
    #     lo_s=event.longitude,
    #     la_r=rec.latitude,
    #     lo_r=rec.longitude,
    #     name=event.name,
    # )
    # plt.savefig(
    #     pjoin(save_folder, f"Location_{event.name}.svg",), dpi=600,
    # )
    # plt.close()

    print(event.name)
    event_nr += 1
    assert event.name == i, "Dictionary and events do not iterate correct"
    # if event.name == "S0235b":
    #     # continue
    #     depths = [38]  # , 32, 59]
    #     # depths = [29]
    # elif event.name == "S0173a":
    #     continue
    #     depths = [14, 38, 65]
    #     # depths = [14]
    # else:
    #     continue
    #     depths = [29, 41]

    if event.name == "S0183a":
        event.distance = 44.5
        event.baz = 73
        event.az = 253
        event.latitude = 15.09
        event.longitude = 179.59
    elif event.name == "S0325a":
        event.distance = 38.4
        event.baz = 0.0
        event.az = 180.0
        from geographiclib.geodesic import Geodesic

        radius = 3389.5
        flattening = 0.0

        dict = Geodesic(a=radius, f=flattening).ArcDirect(
            lat1=rec.latitude, lon1=rec.longitude, azi1=event.baz, a12=event.distance, outmask=1929
        )
        event.latitude = dict["lat2"]
        event.longitude = dict["lon2"]

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
        extra_phases = ["PP", "sP", "pP", "sS", "SS"]

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
        for baz in bazs:
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
            #     plot_extra_phases=extra_phases,
            #     color_plot="blue",
            #     Ylims=ylims,
            #     Parallel=True,
            # )

            # SS_MTI.Inversion.Direct(
            #     fwd=fwd,
            #     misfit=misfit,
            #     event=event,
            #     rec=rec,
            #     phases=phases,
            #     components=components,
            #     phase_corrs=phase_corrs,
            #     t_pre=t_pre,
            #     t_post=t_post,
            #     depths=depths,
            #     tstars=tstars,
            #     fmin=fmin,
            #     fmax=fmax,
            #     zerophase=zerophase,
            #     output_folder=output_folder,
            #     plot=True,
            #     plot_extra_phases=extra_phases,
            #     color_plot="red",
            #     Ylims=ylims,
            #     Parallel=True,
            # )

            """ Post-processing """

            """ (waveform plotting post inversion from generated files)"""
            _PostProcessing.post_waveform_plotting(
                h5_file_folder=output_folder,
                method="GS",
                misfit_name=misfit.name,
                misfit_weight_len=misfit.start_weight_len,
                fwd=fwd,
                event=event,
                rec=rec,
                phases=phases,
                components=components,
                t_pre=t_pre,
                t_post=t_post,
                depths=depths,
                phase_corrs=phase_corrs,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
                tstars=tstars,
                plot_extra_phases=extra_phases,
                Ylims=ylims,
            )

            _PostProcessing.post_waveform_plotting(
                h5_file_folder=output_folder,
                method="Direct",
                misfit_name=misfit.name,
                misfit_weight_len=misfit.start_weight_len,
                fwd=fwd,
                event=event,
                rec=rec,
                phases=phases,
                components=components,
                t_pre=t_pre,
                t_post=t_post,
                depths=depths,
                phase_corrs=phase_corrs,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
                tstars=tstars,
                plot_extra_phases=extra_phases,
                Ylims=ylims,
            )

            # _PostProcessing.post_waveform_plotting_COMBINED(
            #     h5_file_folder=output_folder,
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

            # # """ (misfit vs depth analysis)"""
            DOF = sum([int((x + y) / v["dt"]) for x, y in zip(v["t_pre"], v["t_post"])])
            Moho_d = 24
            # # depths = np.arange(5, 89, 3)
            # fig = _PostProcessing.plot_misfit_vs_depth(
            #     baz=event.baz,
            #     save_paths=[output_folder],
            #     event_name=event.name,
            #     DOF=DOF,
            #     depths=depths,
            #     misfit_name=misfit.name,
            #     veloc_model=fwd.veloc_name,
            #     true_depth=None,
            #     Moho=Moho_d,
            #     fmin=fmin,
            #     fmax=fmax,
            #     amount_of_phases=len(v["phases"]),
            # )
            # plt.tight_layout()
            # plt.savefig(
            #     pjoin(
            #         save_folder,
            #         f"Misfit_vs_Depth_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}_{event.baz}.pdf",
            #     ),
            #     dpi=600,
            # )
            # plt.close()
            # print(f"saved: {save_folder}")

            """ (best MT vs depth phase arrivals) """
            # depths_phases = np.arange(5, 90, 3)  # np.array([23, 26, 29])  #
            # t_pre = [5, 5, 5, 5, 5]
            # t_post = [20, 30, 30, 20, 30]
            # phases = [phases[0], phases[1], phases[2], phases[3], phases[4]]
            # components = [components[0], components[1], components[2], components[3], components[4]]
            # phase_corrs = [
            #     phase_corrs[0],
            #     phase_corrs[1],
            #     phase_corrs[2],
            #     phase_corrs[3],
            #     phase_corrs[4],
            # ]
            # tstars = [tstars[0], tstars[1], tstars[2], tstars[3], tstars[4]]
            # # tstars = [tstar_P, tstar_S]
            # if event.name == "S0173a":
            #     start_depth_range = [32, 77]
            #     end_depth_range = [41, 77]
            # elif event.name == "S0235b":
            #     start_depth_range = [32, 56]
            #     end_depth_range = [32, 65]

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
            #     pref_depth_start=start_depth_range,
            #     pref_depth_end=end_depth_range,
            # )
            # # plt.tight_layout()
            # plt.savefig(
            #     pjoin(
            #         save_folder,
            #         f"Full_PhaseTracking_{event.name}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.svg",
            #     ),
            #     dpi=600,
            # )
            # plt.close()

            # depths_phases = np.arange(5, 90, 3)
            # t_pre = [5, 5]
            # t_post = [30, 30]
            # phases = [phases[0], phases[1]]
            # components = [components[0], components[1]]
            # phase_corrs = [
            #     phase_corrs[0],
            #     phase_corrs[1],
            # ]
            # tstars = [tstars[0], tstars[1]]
            # # tstars = [tstar_P, tstar_S]
            # if event.name == "S0173a":
            #     start_depth_range = [32, 77]
            #     end_depth_range = [41, 77]
            # elif event.name == "S0235b":
            #     start_depth_range = [32, 56]
            #     end_depth_range = [32, 65]

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
            #     pref_depth_start=start_depth_range,
            #     pref_depth_end=end_depth_range,
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
            # depths = np.arange(5, 90, 3)
            # # depts = [14, 29, 32, 35, 41, 38, 65]
            # if event.name == "S0235b":
            #     continue
            #     depts = [29, 32, 35, 38, 56, 59, 62, 65]
            # elif event.name == "S0173a":
            #     depts = [29, 32, 35, 41, 38]
            # fig_sdr = _PostProcessing.Source_Uncertainty(
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

            """ Polarization plot """
            # _PostProcessing.polarization(
            #     h5_file_folder=output_folder,
            #     event=event,
            #     method="GS",
            #     misfit_name=misfit.name,
            #     fwd=fwd,
            #     phases=phases,
            #     components=components,
            #     depths=np.arange(5, 90, 3),
            #     DOF=DOF,
            #     fmin=fmin,
            #     fmax=fmax,
            # )

