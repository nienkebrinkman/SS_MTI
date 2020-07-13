""" 
This script will show source vs receiver affects of certain reflection phases 
"""
__author__ = "Nienke Brinkman"

from obspy.taup import TauPyModel
import obspy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os
import instaseis

import SS_MTI
from EventInterface import EventObj
import SS_MTI.DataGetter as DG
import SS_MTI.PostProcessing as PP


""" Define parameters """
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 10.99032013
lon_src = 170
name = "Test_Event"

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

epi, az, baz = EventObj.Get_location(lat_src, lon_src, lat_rec, lon_rec)

fmin = 1.0 / 10.0
fmax = 1.0 / 3.0
zerophase = False


""" Define different velocity models"""
db_name_1 = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
npz_file_name_1 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK_BKE.npz"

db_name_2 = "/mnt/marshost/instaseis2/databases/TAYAK_shallow"
npz_file_name_2 = "/home/nienke/Documents/Research/Data/npz_files/TAYAK.npz"

db_name_3 = "/mnt/marshost/instaseis/databases/blindtestmodels_1s/DWThot_1s"
npz_file_name_3 = "/home/nienke/Documents/Research/Data/npz_files/DWThot.npz"

db_names = [db_name_1, db_name_2]
npz_file_names = [npz_file_name_1, npz_file_name_2]


mnt_folder = "/mnt/marshost/"

# DG.unmnt_remote_folder(mnt_folder=mnt_folder)
DG.mnt_remote_folder(
    host_ip="marshost.ethz.ch", host_usr="sysop", remote_folder="/data/", mnt_folder=mnt_folder,
)

""" Define the depths you want to loop over """
depths = np.arange(5, 90, 10)
# depths = [5, 30]

for npz_file, db_path in zip(npz_file_names, db_names):

    model = TauPyModel(npz_file)
    db = instaseis.open_db(db_path)

    model_name = npz_file.split("/")[-1].strip(".npz")
    folder = "/home/nienke/Documents/Research/Data/MTI/Phase_tracking"
    save_path = os.path.join(folder, model_name)
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    veloc_fig = PP.Plot_veloc_models(model)
    veloc_fig.savefig(os.path.join(save_path, "Velocity_model.pdf"))

    """ Get all the reflection phases """
    Moho = str(int(model.model.moho_depth))

    Underside_refl_src = []
    Conversion_src = []
    Conversion_rec = []
    for i, values in enumerate(model.model.s_mod.critical_depths):
        if values[0] < 1.0 or values[0] > 100.0:
            pass
        else:
            interface = str(int(values[0]))
            for down_phase in ["p^", "s^"]:
                for up_phase in ["P", "S"]:
                    Underside_refl_src.append(down_phase + interface + up_phase)
            Conversion_src.append("S" + interface + "P")
            Conversion_src.append("P" + interface + "S")
            # Conversion_rec.append("P" + interface + "p")
            Conversion_rec.append("P" + interface + "s")

    Direct_phases = ["P", "S"]
    Depth_phases = ["pP", "sP", "sS"]
    Double_phases = ["PP", "PPP", "SSS"]
    phases = (
        Direct_phases
        + Depth_phases
        + Double_phases
        + Conversion_src
        + Conversion_rec
        + Underside_refl_src
    )
    phase_colors = (
        ["grey"] * len(Direct_phases)
        + ["grey"] * len(Depth_phases)
        + ["grey"] * len(Double_phases)
        + ["red"] * len(Conversion_src)
        + ["blue"] * len(Conversion_rec)
        + ["green"] * len(Underside_refl_src)
    )
    phase_labels = {
        "grey": "Direct/Double/Depth phase",
        "red": "source side",
        "blue": "receiver side",
        "green": "underside reflection",
    }

    """ Vizualize the reflection phases """
    # if not os.path.exists(os.path.join(save_path, "ray_paths")):
    #     os.mkdir(os.path.join(save_path, "ray_paths"))
    # for depth in depths:
    #     arrivals = model.get_ray_paths(
    #         source_depth_in_km=depth, distance_in_degree=epi, phase_list=Conversion_src
    #     )

    #     ax = arrivals.plot_rays(plot_type="cartesian", show=False, legend=True)
    #     plt.savefig(os.path.join(save_path, "ray_paths", f"src_d_{depth}"))
    #     plt.close()

    #     arrivals = model.get_ray_paths(
    #         source_depth_in_km=depth, distance_in_degree=epi, phase_list=Conversion_rec
    #     )

    #     ax = arrivals.plot_rays(plot_type="cartesian", show=False, legend=True)
    #     plt.savefig(os.path.join(save_path, "ray_paths", f"rec_d_{depth}"))

    """ Define forward solver """

    fwd = SS_MTI.Forward.Instaseis(
        instaseis_db=db,
        taup_model=npz_file,
        rec_lat=lat_rec,
        rec_lon=lon_rec,
        or_time=or_time,
        dt=dt,
        start_cut=0.0,
        end_cut=800.0,
    )

    figP = None
    axP = None
    figS = None
    axS = None
    Yticks = np.arange(len(depths))

    normalize = True

    for Ytick, depth in zip(Yticks, depths):
        """ GENERATE GREEN'S FUNCTION AT SPECIFIC DEPTH """
        syn_tts = []
        syn_GFs = []
        for phase in phases:
            print(phase)
            syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=epi)
            syn_tts.append(syn_tt)

        st_syn = obspy.Stream()
        for i, comp in enumerate(list(components)):
            """ Compute GreensFunctions at depth"""
            syn_GF = fwd.get_greens_functions(
                comp=comp,
                depth=depth,
                distance=epi,
                lat_src=lat_src,
                lon_src=lon_src,
                tstar=None,
                LQT=False,
                inc=None,
                baz=baz,
                M0=M0,
                filter=True,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
            )
            syn_GFs.append(syn_GF)

            """ Generate synthetic waveform"""
            tr_syn = fwd.generate_synthetic_data(
                st_GF=syn_GFs[i], focal_mech=focal_mech, M0=M0, slice=False,
            )
            tr_syn.stats.channel = f"BX{comp}"
            st_syn += tr_syn

        """ Normalize the stream"""
        if normalize:
            st_syn.normalize(global_max=True)

        extra_phases = phases[2:]
        extra_arrs = syn_tts[2:]

        """ Phase vs Depth """

        ## P - Phase
        figP, axP = PP.Plot_trace_vs_depth(
            stream=st_syn[0:2],
            depth=depth,
            total_depths=len(depths),
            Ytick=Ytick,
            phase=Direct_phases[0],
            phase_arr=syn_tts[0],
            t_pre=10.0,
            t_post=35.0,
            fig=figP,
            ax=axP,
            extra_phases=extra_phases,
            extra_arrs=extra_arrs,
            phase_colors=phase_colors,
            phase_labels=phase_labels,
        )

        ## S - Phase
        figS, axS = PP.Plot_trace_vs_depth(
            stream=st_syn,
            depth=depth,
            total_depths=len(depths),
            Ytick=Ytick * 2,
            phase=Direct_phases[1],
            phase_arr=syn_tts[1],
            t_pre=10.0,
            t_post=50.0,
            fig=figS,
            ax=axS,
            extra_phases=extra_phases,
            extra_arrs=extra_arrs,
            phase_colors=phase_colors,
            phase_labels=phase_labels,
        )

        """ Phases at one depth """
        fig1, ax1 = PP.Plot_phases_vs_comp(
            stream=st_syn,
            phase_cuts=Direct_phases,
            phase_arrs=syn_tts[0 : len(Direct_phases)],
            t_pre=10.0,
            t_post=50.0,
            extra_phases=extra_phases,
            extra_arrs=extra_arrs,
            phase_colors=phase_colors,
            phase_labels=phase_labels,
        )
        fig1.text(
            0.55,
            0.9,
            f"Veloc model: {model_name},Depth:{depth}",
            ha="center",
            va="bottom",
            size="x-large",
            color="blue",
        )
        fig1.savefig(os.path.join(save_path, f"d_{depth}.png"))

    axS[0].yaxis.set_ticks(Yticks * 2)
    axS[0].set_yticklabels(depths)
    figS.text(
        0.5,
        0.95,
        f"Velocity model: {model_name}",
        ha="center",
        va="bottom",
        size="x-large",
        color="blue",
    )
    figS.savefig(os.path.join(save_path, f"{Direct_phases[1]}_Phase_model_{model_name}.pdf"))

    axP[0].yaxis.set_ticks(Yticks)
    axP[0].set_yticklabels(depths)
    figP.text(
        0.5,
        0.95,
        f"Velocity model: {model_name}",
        ha="center",
        va="bottom",
        size="x-large",
        color="blue",
    )
    figP.savefig(os.path.join(save_path, f"{Direct_phases[0]}_Phase_model_{model_name}.pdf"))

DG.unmnt_remote_folder(mnt_folder=mnt_folder)
