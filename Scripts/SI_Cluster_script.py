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
import argparse
import toml
import mpi4py.MPI
from os import makedirs

import SS_MTI
import EventInterface
from SS_MTI import PostProcessing as _PostProcessing


def define_arguments():
    helptext = "Determine focal mechanisms of Marsquake"
    parser = argparse.ArgumentParser(description=helptext)

    helptext = "Input toml file"
    parser.add_argument("input_file", help=helptext)
    return parser.parse_args()


if __name__ == "__main__":
    Parallel = True
    if Parallel:
        print("Your inversion will be run in parallel")

    # input_file = "/home/nienke/Documents/Research/SS_MTI/Input/TAYAK_BKE_tstar_update.toml"
    args = define_arguments()
    print(f"Inversion based on input file: {args.input_file}")
    event_input = toml.load(args.input_file, _dict=dict)
    save_folder = pjoin(
        "/home/nienke/Data_2020/Test_2021/", args.input_file.split("/")[-1].strip(".toml")
    )
    if not exist(save_folder):
        makedirs(save_folder)

    # save_folder = "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_2/Test/"

    path = "/home/nienke/Data_2020/catalog"
    # path = "/home/nienke/Documents/Research/Data/MTI/catalog"
    path_to_inventory = pjoin(path, "inventory.xml")
    path_to_catalog = pjoin(path, "catalog.xml")

    """ Read the inventory and catalog file (the once that contain info about the marsquakes) """
    inv = None  # SS_MTI.DataGetter.read_inv(inv_path=path_to_inventory)  # Inventory file
    cat = SS_MTI.DataGetter.read_cat(cat_path=path_to_catalog)  # Catalog file

    """ Get the data into a list of obspy.Event objects """
    events = SS_MTI.DataGetter.read_events_from_cat(
        event_params=event_input,
        cat=cat,
        inv=inv,
        local_folder=pjoin(save_folder, "event.mseed"),
        host_name=None,
        user_name=None,
        remote_folder=None,
        save_file_name=pjoin(save_folder, "event.mseed"),
    )

    """ Specify receiver """
    lat_rec = 4.5  # 02384
    lon_rec = 135.623447
    rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

    """ """
    depths = np.arange(5, 90, 3)
    # depths = np.arange(29, 50, 3)
    # depths = [59]

    strikes = np.arange(0, 360, 15)
    dips = np.arange(0, 91, 10)
    rakes = np.arange(-180, 180, 15)

    # strikes = [15.0116557194]  # [132.395557582]
    # dips = [59.551091053]  # [51.9591191063]
    # rakes = [-45.6275510954]  # [-139.94976385]

    # bazs = np.arange(0, 360, 20)

    """ Define different velocity models"""
    db_name_1 = "/opt/databases/TAYAK_15s_BKE"
    npz_file_name_1 = "/home/nienke/Data_2020/npz_files/TAYAK_BKE.npz"

    db_name_2 = "/opt/databases/TAYAK_shallow"
    npz_file_name_2 = "/home/nienke/Data_2020/npz_files/TAYAK.npz"

    # db_name_1 = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
    # npz_file_name_1 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz"

    # db_name_2 = "/mnt/marshost/instaseis2/databases/TAYAK_shallow"
    # npz_file_name_2 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"

    # db_name_3 = "/mnt/marshost/instaseis2/databases/TAYAK_1s_30km"
    # npz_file_name_3 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_30km.npz"

    db_names = [db_name_1]  # , db_name_3, db_name_4, db_name_5]
    npz_file_names = [npz_file_name_1]

    """ Loop over events to invert for: """
    event_nr = 0
    for i, v in event_input.items():
        event = events[event_nr]
        print(event.name)
        event_nr += 1
        assert event.name == i, "Dictionary and events do not iterate correct"
        # if event.name == "S0173a" or event.:
        #     pass
        # else:
        #     continue

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
                lat1=rec.latitude,
                lon1=rec.longitude,
                azi1=event.baz,
                a12=event.distance,
                outmask=1929,
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
            # mnt_folder = "/mnt/marshost/"
            # if not lsdir(mnt_folder):
            #     print(f"{mnt_folder} is still empty, mounting now...")
            #     SS_MTI.DataGetter.mnt_remote_folder(
            #         host_ip="marshost.ethz.ch",
            #         host_usr="sysop",
            #         remote_folder="/data/",
            #         mnt_folder=mnt_folder,
            #     )

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
                misfit = SS_MTI.Misfit.L2(
                    weights=weights, start_weight_len=start_weight_len, dt=dt
                )
            elif misfit_method == "CC":
                misfit = SS_MTI.Misfit.CC(shift_samples=128)
            elif misfit_method == "POL":
                misfit = SS_MTI.Misfit.Pol(
                    components=components,
                    start_weight_len=start_weight_len,
                    weights=weights,
                    dt=dt,
                )
            else:
                raise ValueError("misfit can be L2, CC or POL in [MISFIT] of .toml file")

            """ Start inversion """
            # for baz in bazs:
            #     event.baz = baz
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
                Parallel=Parallel,
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
                Parallel=Parallel,
            )

            """ Post-processing """
            # if Parallel:
            #     mpi4py.MPI.COMM_WORLD.Barrier()
            #     rank = mpi4py.MPI.COMM_WORLD.Get_rank()
            #     if not rank == 0:
            #         print(
            #             "rank {rank} will go to next simulation and does not doe post processing"
            #         )
            #         continue

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

            """ (misfit vs depth analysis)"""
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
            # t_post = [30, 30]
            # phases = [phases[0], phases[1]]
            # components = [components[0], components[1]]
            # phase_corrs = [phase_corrs[0], phase_corrs[1]]
            # tstars = [tstars[0], tstars[1]]
            # # tstars = [tstar_P, tstar_S]
            # start_depth_range = 29  # 53  #
            # end_depth_range = 41  # 69  #
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
            # fig_sdr = _PostProcessing.Source_Uncertainty(
            #     h5_file_folder=output_folder,
            #     event_name=event.name,
            #     method="GS",
            #     misfit_name=misfit.name,
            #     fwd=fwd,
            #     phases=phases,
            #     components=components,
            #     depths=np.arange(53, 68, 3),
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

