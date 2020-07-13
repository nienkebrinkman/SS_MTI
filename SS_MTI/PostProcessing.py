import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import obspy


def Plot_veloc_models(Taup_model, depth_event=None, depth_syn=None):
    depth = np.array([])
    Vp = np.array([])
    Vs = np.array([])
    dens = np.array([])

    for i, values in enumerate(Taup_model.model.s_mod.v_mod.layers):
        depth = np.append(depth, values[0])
        depth = np.append(depth, values[1])
        Vp = np.append(Vp, values[2])
        Vp = np.append(Vp, values[3])
        Vs = np.append(Vs, values[4])
        Vs = np.append(Vs, values[5])
        dens = np.append(dens, values[6])
        dens = np.append(dens, values[7])

    fig, ax = plt.subplots(1, 3, sharey="all", sharex="all", figsize=(8, 6))
    ax[0].plot(Vp, depth)
    if depth_event is not None:
        int_vp = interpolate.interp1d(depth, Vp)
        event_vp = int_vp(depth_event)
        ax[0].plot(event_vp, depth_event, "g*", markersize=15, label="Event Depth")
        if depth_syn is not None:
            for i in range(len(depth_syn)):
                event_vp = int_vp(depth_syn[i])
                if i == 0:
                    ax[0].plot(
                        event_vp, depth_syn[i], "r*", markersize=15, label="Synthetic Depth"
                    )
                else:
                    ax[0].plot(event_vp, depth_syn[i], "r*", markersize=15, label="_hidden")
    ax[0].set_title("VP", color="b", fontsize=20)
    ax[0].set_ylabel("Depth [km]", fontsize=20)
    ax[0].tick_params(axis="x", labelsize=18)
    ax[0].tick_params(axis="y", labelsize=18)
    # ax[0].ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
    ax[0].grid(True)
    # ax[0].set_ylim([500,0])
    ax[0].set_xlim([0, 8])

    ax[1].plot(Vs, depth, label="Shallow")
    if depth_event is not None:
        int_vs = interpolate.interp1d(depth, Vs)
        event_vs = int_vs(depth_event)
        ax[1].plot(event_vs, depth_event, "g*", markersize=15, label="Event Depth")
        if depth_syn is not None:
            for i in range(len(depth_syn)):
                event_vs = int_vs(depth_syn[i])
                if i == 0:
                    ax[1].plot(
                        event_vs, depth_syn[i], "r*", markersize=15, label="Synthetic Depth"
                    )
                else:
                    ax[1].plot(event_vs, depth_syn[i], "r*", markersize=15, label="_hidden")
    # ax[1].legend()
    ax[1].set_title("VS", color="b", fontsize=20)
    ax[1].tick_params(axis="x", labelsize=18)
    ax[1].tick_params(axis="y", labelsize=18)
    # ax[1].ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
    ax[1].grid(True)
    # ax[0].set_ylim([0,100])

    ax[2].plot(dens, depth)
    if depth_event is not None:
        int_dens = interpolate.interp1d(depth, dens)
        event_dens = int_dens(depth_event)
        ax[2].plot(event_dens, depth_event, "g*", markersize=15, label="Event Depth")
        if depth_syn is not None:
            for i in range(len(depth_syn)):
                event_dens = int_dens(depth_syn[i])
                if i == 0:
                    ax[2].plot(
                        event_dens, depth_syn[i], "r*", markersize=15, label="Synthetic Depth"
                    )
                else:
                    ax[2].plot(event_dens, depth_syn[i], "r*", markersize=15, label="_hidden")
        ax[2].legend()
    ax[2].set_title("Density", color="b", fontsize=20)
    ax[2].tick_params(axis="x", labelsize=18)
    ax[2].tick_params(axis="y", labelsize=18)
    # ax[2].ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
    ax[2].grid(True)
    ax[2].set_ylim([0, 100])
    ax[0].set_ylim(ax[0].get_ylim()[::-1])
    return fig


def Plot_trace_vs_depth(
    stream: obspy.Stream,
    depth: float,
    total_depths: int,
    Ytick: float,
    phase: str,
    phase_arr: float,
    t_pre: float = 10.0,
    t_post: float = 50.0,
    fig: plt.figure = None,
    ax: plt.axes = None,
    extra_phases: [str] = None,
    extra_arrs: [float] = None,
    phase_colors: [str] = None,
    phase_labels: dict = None,
):
    if fig is None and ax is None:
        fig, ax = plt.subplots(
            nrows=1,
            ncols=len(stream),
            figsize=(5 * len(stream), 2 * total_depths),
            sharex="col",
            sharey="all",
        )
    st = stream.copy()
    st.trim(
        starttime=st[0].stats.starttime + phase_arr - t_pre,
        endtime=st[0].stats.starttime + phase_arr + t_post,
    )

    global_max = max([tr.data.max() for tr in st])
    global_min = min([tr.data.min() for tr in st])
    y = global_max * 0.9 + Ytick

    ymin = global_min + Ytick
    ymax = global_max + Ytick
    for i in range(len(stream)):
        ax[i].plot(
            st[i].times() - t_pre, st[i].data + Ytick, "k",
        )

        ax[i].plot(
            [0, 0], [ymin, ymax], "grey",
        )
        ax[i].text(0, y, phase, verticalalignment="center", color="grey", fontsize=6)
        if extra_phases is not None:
            for k in range(len(extra_phases)):
                if extra_arrs[k] is None:
                    continue
                phase_t = extra_arrs[k] - phase_arr
                if phase_colors is None:
                    c = "grey"
                else:
                    c = phase_colors[k]
                ax[i].plot(
                    [phase_t, phase_t], [ymin, ymax], c,
                )
                ax[i].text(
                    phase_t + 0.1,
                    y,
                    extra_phases[k],
                    verticalalignment="center",
                    color=c,
                    fontsize=6,
                    rotation=90,
                )
        ax[i].set_xlim(-t_pre, t_post)
        ax[i].set_title(f"{phase}-Phase channel:{st[i].stats.channel}")
    if phase_colors is not None:
        unique_colors = list(set(phase_colors))
        # unique_list = [mpatches.Patch(color=c, label=phase_labels[c]) for c in phase_labels]
        unique_list = [
            Line2D([0], [0], color=c, linewidth=3, label=phase_labels[c]) for c in phase_labels
        ]
    ax[0].legend(
        handles=unique_list, prop={"size": 6}, loc="upper left", bbox_to_anchor=(0.0, 1.07)
    )
    # fig.legend(handles=unique_list, prop={"size": 6}, loc="upper left")
    fig.text(0.04, 0.5, "Source Depth (km)", va="center", rotation="vertical")
    fig.text(0.5, 0.04, "Time after arrival (s)", va="center")
    return fig, ax


def Plot_phases_vs_comp(
    stream: obspy.Stream,
    phase_cuts: [str],
    phase_arrs: [float],
    t_pre: float = 20.0,
    t_post: float = 60.0,
    extra_phases: [str] = None,
    extra_arrs: [float] = None,
    phase_colors: [str] = None,
    phase_labels: dict = None,
):
    """ Plotting function that cuts the stream the phases in phase_cuts"""
    if not len(phase_cuts) == len(phase_arrs):
        raise ValueError("phase_cut and phase_arrs should have same length")
    if extra_phases is not None:
        if not len(extra_phases) == len(extra_arrs):
            raise ValueError("extra_phases and extra_arrs should have same length")
    fig, ax = plt.subplots(
        nrows=len(stream), ncols=len(phase_cuts), figsize=(18, 8), sharex="col", sharey="all"
    )
    for j in range(len(phase_cuts)):

        st = stream.copy()
        st.trim(
            starttime=st[0].stats.starttime + phase_arrs[j] - t_pre,
            endtime=st[0].stats.starttime + phase_arrs[j] + t_post,
        )

        for i in range(len(stream)):
            ax[i, j].plot(
                st[i].times() - t_pre, st[i].data, "k",
            )
            y = ax[i, j].get_ylim()[1] * 0.8
            ax[i, j].axvline(x=0, c="grey")
            ax[i, j].text(
                0, y, phase_cuts[j], verticalalignment="center", color="grey", fontsize=6
            )
            if extra_phases is not None:
                for k in range(len(extra_phases)):
                    if extra_arrs[k] is None:
                        continue

                    if phase_colors is None:
                        c = "grey"
                    else:
                        c = phase_colors[k]
                    phase_t = extra_arrs[k] - phase_arrs[j]
                    ax[i, j].axvline(x=phase_t, c=c)
                    ax[i, j].text(
                        phase_t + 0.1,
                        y,
                        extra_phases[k],
                        verticalalignment="center",
                        color=c,
                        fontsize=6,
                        rotation=90,
                    )
            ax[i, j].set_xlim(-t_pre, t_post)
            if i == 0:
                ax[i, j].set_title(f"{phase_cuts[j]}-phase")

            if j == 0:
                ax[i, j].set_ylabel(st[i].stats.channel)
    ax[0, 0].set_ylim(-1, 1)
    if phase_colors is not None:
        unique_colors = list(set(phase_colors))
        # unique_list = [mpatches.Patch(color=c, label=phase_labels[c]) for c in phase_labels]
        unique_list = [
            Line2D([0], [0], color=c, linewidth=3, label=phase_labels[c]) for c in phase_labels
        ]
    ax[0, 0].legend(
        handles=unique_list, prop={"size": 6}, loc="upper left", bbox_to_anchor=(0.0, 1.4)
    )
    fig.text(0.04, 0.5, "Displacement (m)", va="center", rotation="vertical")
    fig.text(0.5, 0.04, "Time after arrival (s)", va="center")
    return fig, ax
