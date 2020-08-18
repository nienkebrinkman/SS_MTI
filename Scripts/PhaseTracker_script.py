""" 
This script will show source vs receiver affects of certain reflection phases 
"""
__author__ = "Nienke Brinkman"

from obspy.taup import TauPyModel
import obspy
import numpy as np
import matplotlib.pyplot as plt


import os
import instaseis

import SS_MTI
from EventInterface import EventObj
import SS_MTI.DataGetter as DG
import SS_MTI.PostProcessing as PP


""" Define parameters """
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 10.99032013
lon_src = 160.94
name = "Test_Event"

lat_rec = 4.502384
lon_rec = 135.623447
rec = instaseis.Receiver(latitude=lat_rec, longitude=lon_rec)

fig = PP.Plot_event_location(lat_src, lon_src, lat_rec, lon_rec, "test event")
plt.savefig("/home/nienke/Documents/Research/Data/MTI/Phase_tracking/source_location.pdf")
plt.close()

strike = 70
dip = 30
rake =-95
focal_mech = [strike, dip, rake]
M0 = 5.62e13

dt = 0.05

components = "ZRT"
kind = "displacement"
noise = True

epi, az, baz = EventObj.Get_location(lat_src, lon_src, lat_rec, lon_rec)
epi = 25.89
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

db_name_4 = "/mnt/marshost/instaseis/databases/blindtestmodels_1s/DWAK_1s"
npz_file_name_4 = "/home/nienke/Documents/Research/Data/npz_files/DWAK.npz"

db_name_5 = "/mnt/marshost/instaseis/databases/blindtestmodels_1s/MAAK_1s"
npz_file_name_5 = "/home/nienke/Documents/Research/Data/npz_files/MAAK.npz"

db_names = [db_name_1]#, db_name_2, db_name_3, db_name_4, db_name_5]
npz_file_names = [
    npz_file_name_1]#,
#     npz_file_name_2,
#     npz_file_name_3,
#     npz_file_name_4,
#     npz_file_name_5,
# ]

mnt_folder = "/mnt/marshost/"

# DG.unmnt_remote_folder(mnt_folder=mnt_folder)
DG.mnt_remote_folder(
    host_ip="marshost.ethz.ch", host_usr="sysop", remote_folder="/data/", mnt_folder=mnt_folder,
)

""" Define the depths you want to loop over """
depths = np.arange(5, 90, 10)
# depths = [85]

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
    Reflection_phases = []
    for i, values in enumerate(model.model.s_mod.critical_depths):
        if values[0] < 1.0 or values[0] > 100.0:
            continue
        interface = str(int(values[0]))
        if i > 1:
            Reflection_phases.append(
                "P" + interface + "s" + str(int(model.model.s_mod.critical_depths[i - 1][0])) + "p"
            )
        if values[0] > 50.0 and "TAYAK" in model_name:
            continue
        for down_phase in ["p^", "s^"]:
            for up_phase in ["P", "S"]:

                Underside_refl_src.append(down_phase + interface + up_phase)

        Conversion_src.append("S" + interface + "P")
        Conversion_src.append("P" + interface + "S")
        Conversion_rec.append("P" + interface + "p")
        Conversion_rec.append("P" + interface + "s")

    Direct_phases = ["P"]
    Depth_phases = ["pP", "sP", "sS"]
    Double_phases = ["PP", "PPP", "SSS"]
    extra_phases = (
        Depth_phases
        + Double_phases
        + Conversion_src
        + Conversion_rec
        + Underside_refl_src
        + Reflection_phases
    )
    extra_phase_colors = (
        ["grey"] * len(Depth_phases)
        + ["grey"] * len(Double_phases)
        + ["red"] * len(Conversion_src)
        + ["blue"] * len(Conversion_rec)
        + ["green"] * len(Underside_refl_src)
        + ["purple"] * len(Reflection_phases)
    )
    phase_labels = {
        "grey": "Direct/Double/Depth phase",
        "red": "source side",
        "blue": "receiver side",
        "green": "underside reflection",
    }

    """ Vizualize the reflection phases """
    if not os.path.exists(os.path.join(save_path, "ray_paths")):
        os.mkdir(os.path.join(save_path, "ray_paths"))
        for depth in depths:
            for phase in extra_phases:
                arrivals = model.get_ray_paths(
                    source_depth_in_km=depth, distance_in_degree=epi, phase_list=[phase]
                )
                if not arrivals:
                    continue
                ax = arrivals.plot_rays(plot_type="cartesian", show=False, legend=True)
                plt.savefig(os.path.join(save_path, "ray_paths", f"d_{depth}_{phase}"))
                plt.close()

    """ Define forward solver """

    fwd = SS_MTI.Forward.Instaseis(
        instaseis_db=db, taup_model=npz_file, or_time=or_time, dt=dt, start_cut=0.0, end_cut=800.0,
    )

    fig = [None] * len(Direct_phases)
    ax = [None] * len(Direct_phases)
    Yticks = np.arange(len(depths)) * 1.8

    normalize = True
    vlines = False

    extra_arrs = [[] for _ in range(len(depths))]

    for idepth, depth in enumerate(depths):
        print(f"Depth {depth}")

        """ Get phase arrivals (not yet the extra phases) """
        phase_arrs = []
        for phase in Direct_phases:
            arr = fwd.get_phase_tt(phase=phase, depth=depth, distance=epi)
            phase_arrs.append(arr)

        """ Generate Green's function at the current depth """
        syn_GFs = []
        st_syn_full = obspy.Stream()
        for i, comp in enumerate(list(components)):
            """ Compute GreensFunctions at depth"""
            syn_GF = fwd.get_greens_functions(
                comp=comp,
                depth=depth,
                distance=epi,
                rec=rec,
                lat_src=lat_src,
                lon_src=lon_src,
                tstar=None,
                LQT=False,
                inc=None,
                baz=baz,
                M0=M0,
                filter=False,
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
            st_syn_full += tr_syn

        for iphase, phase in enumerate(Direct_phases):
            for extraphase in extra_phases:
                arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=epi)
                if arr:
                    extra_arrs[idepth].append(arr - phase_arrs[iphase])
                else:
                    extra_arrs[idepth].append(arr)
            st_syn = st_syn_full.copy()

            t_pre = 5.0
            t_post = 35.0

            """ Trim the stream around the phase """
            st_syn.trim(
                starttime=st_syn[0].stats.starttime + phase_arrs[iphase] - t_pre,
                endtime=st_syn[0].stats.starttime + phase_arrs[iphase] + t_post,
            )

            """ Normalize the stream"""
            if normalize:
                st_syn.normalize(global_max=True)

            """ Phase vs Depth """
            if phase == "P":
                stream = st_syn[0:2]
            else:
                stream = st_syn

            if vlines:
                fig[iphase], ax[iphase] = PP.Plot_trace_vs_depth_copy(
                    stream=stream,
                    depth=depth,
                    total_depths=len(depths),
                    Ytick=Yticks[idepth],
                    phase=phase,
                    phase_arr=phase_arrs[iphase],
                    t_pre=t_pre,
                    t_post=t_post,
                    fig=fig[iphase],
                    ax=ax[iphase],
                    extra_phases=extra_phases,
                    extra_arrs=extra_arrs[idepth],
                    phase_colors=extra_phase_colors,
                    phase_labels=phase_labels,
                )
            else:
                fig[iphase], ax[iphase] = PP.Plot_trace_vs_depth(
                    stream=stream,
                    phase=phase,
                    total_depths=len(depths),
                    Ytick=Yticks[idepth],
                    t_pre=t_pre,
                    t_post=t_post,
                    fig=fig[iphase],
                    ax=ax[iphase],
                )
    delta = Yticks[1] - Yticks[0]
    for iphase, phase in enumerate(Direct_phases):
        if vlines:
            pass
        else:
            for j in range(len(list(components))):
                if phase == "P" and j == 2:
                    continue
                for k in range(len(extra_phases)):
                    x = np.asarray([arr[k] for arr in extra_arrs], dtype=np.float)
                    y = np.asarray(Yticks)
                    y = y[~np.isnan(x)]
                    x = x[~np.isnan(x)]
                    if x.size == 0:
                        continue
                    rotn = np.degrees(np.arctan(y[-1:] - y[-2:-1], x[-1:] - x[-2:-1]))
                    if rotn.size == 0:
                        trans_angle = 90
                        ax[iphase][j].plot(
                            [x[0], x[0]], [y[0] - 0.8, y[0] + 0.8], extra_phase_colors[k],
                        )
                    else:
                        l2 = np.array((x[-1], y[-1]))
                        rotation = rotn[-1]
                        trans_angle = plt.gca().transData.transform_angles(
                            np.array((rotation,)), l2.reshape((1, 2))
                        )[0]
                        ax[iphase][j].plot(x, y, "-", c=extra_phase_colors[k])

                    ax[iphase][j].text(
                        x[-1],
                        y[-1],
                        extra_phases[k],
                        verticalalignment="center",
                        color=extra_phase_colors[k],
                        fontsize=6,
                        rotation=trans_angle,
                    )

        ax[iphase][0].yaxis.set_ticks(Yticks)
        ax[iphase][0].set_yticklabels(depths)
        ax[iphase][0].set_ylim(Yticks[0] - delta, Yticks[-1] + delta)
        fig[iphase].text(
            0.5,
            0.95,
            f"Velocity model: {model_name}",
            ha="center",
            va="bottom",
            size="x-large",
            color="blue",
        )
        fig[iphase].savefig(os.path.join(save_path, f"{phase}_Phase_model_{model_name}.pdf"))
        plt.close(fig[iphase])

        a = 1

DG.unmnt_remote_folder(mnt_folder=mnt_folder)
