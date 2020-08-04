import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import obspy
from PIL import Image
from matplotlib.lines import Line2D
import os
from matplotlib.patches import Circle
from obspy.imaging.beachball import beach
import matplotlib.image as mpimg
import io
from sklearn import preprocessing
import pandas as pd
import matplotlib.cm as cm

pyproj_datadir = os.environ["PROJ_LIB"]

from mpl_toolkits.basemap import Basemap
import re


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


def Plot_trace_vs_depth_copy(
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
                phase_t = extra_arrs[k]
                if phase_colors is None:
                    y = global_max * 0.9 + Ytick
                    c = "grey"
                else:
                    y = global_max * 0.9 + Ytick
                    ind = re.findall(r"\d+", extra_phases[k])
                    if ind:
                        if len(ind) == 2:
                            if int(ind[0]) < depth:
                                c = "blue"
                                y = global_min * 0.4 + Ytick
                            else:
                                c = "red"
                                y = global_min * 0.4 + Ytick
                        else:
                            c = phase_colors[k]
                    else:
                        c = phase_colors[k]
                        y = global_max * 0.9 + Ytick

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


def Plot_trace_vs_depth(
    stream: obspy.Stream,
    phase: str,
    total_depths: int,
    Ytick: float,
    t_pre: float = 10.0,
    t_post: float = 50.0,
    fig: plt.figure = None,
    ax: plt.axes = None,
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
    global_max = max([tr.data.max() for tr in st])
    global_min = min([tr.data.min() for tr in st])
    y = global_max * 0.9 + Ytick
    ymin = global_min + Ytick
    ymax = global_max + Ytick
    for i in range(len(stream)):
        ax[i].plot(
            st[i].times() - t_pre, st[i].data + Ytick, "k",
        )
        ax[i].set_xlim(-t_pre, t_post)
        ax[i].set_title(f"{phase}-Phase channel:{st[i].stats.channel}")
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


def Plot_event_location(
    la_s: float, lo_s: float, la_r: float, lo_r: float, name: str = "test event"
):
    # la_s = event.latitude
    # lo_s = event.longitude

    mars_dir = "/home/nienke/Documents/Research/Data/mars_pictures/Mars_lightgray.jpg"

    fig = plt.figure(figsize=(10, 8))

    # m = Basemap(projection='moll', lon_0=round(0.0))
    m = Basemap(
        projection="merc", llcrnrlat=-80, urcrnrlat=80, llcrnrlon=0, urcrnrlon=200, resolution="c"
    )

    # draw parallels and meridians.
    par = np.arange(-90, 90, 30)
    label_par = np.full(len(par), True, dtype=bool)
    meridians = np.arange(-180, 180, 30)
    label_meri = np.full(len(meridians), True, dtype=bool)

    m.drawmeridians(np.arange(-180, 180, 30), labels=label_meri)
    m.drawparallels(np.arange(-90, 90, 30), label=label_par)

    m.warpimage(mars_dir)
    mstatlon, mstatlat = m(lo_r, la_r)
    m.plot(mstatlon, mstatlat, "k^", markersize=20, label="InSight")

    EQlonA, EQlatA = m(lo_s, la_s)
    #     EQlonB, EQlatB = m(lo_sB, la_sB)
    #     EQlonC, EQlatC = m(lo_sC, la_sC)
    m.plot(EQlonA, EQlatA, "r*", markersize=20, zorder=10, label=name)
    # m.plot(EQlonB, EQlatB, 'g*', markersize=20, zorder=10, label = event_B.name)
    #     m.plot(EQlonC, EQlatC, 'b*', markersize=20, zorder=10, label=event_C.name)
    plt.legend(fontsize=20)
    plt.tight_layout()
    # plt.show()
    #     plt.savefig('Location_Event.pdf')
    return fig


def Plot_waveforms(
    fig: plt.figure, ax: plt.axes,
):
    pass


""" Plot beachballs """


def Get_bb_img(MT, color, alpha=1.0):
    ### FULL MOMENT TENSOR
    img = None
    buf = io.BytesIO()
    fig_bb = plt.figure(figsize=(5, 5), dpi=200)
    ax_bb_1 = fig_bb.add_axes([0.0, 0.0, 1.0, 1.0])
    ax_bb_1.set_xticks([])
    ax_bb_1.set_yticks([])
    ax_bb_1.axis("off")

    if np.count_nonzero(MT) < 6 and len(MT) == 6:
        pass
    else:
        b = beach(
            fm=MT, width=990, linewidth=0, facecolor=color, xy=(0, 0), axes=ax_bb_1, alpha=alpha
        )
        ax_bb_1.add_collection(b)
    ax_bb_1.set_xlim((-1, 1))
    ax_bb_1.set_ylim((-1, 1))

    buf.seek(0)
    fig_bb.savefig(buf, format="png", dpi=200)
    buf.seek(0)
    if img is None:
        img = mpimg.imread(buf)
    else:
        img += mpimg.imread(buf)

    plt.close(fig_bb)
    return img, buf


def Plot_Direct_BB(
    MT_Full,
    Eps,
    MT_DC,
    M0_DC,
    MT_CLVD,
    M0_CLVD,
    azimuths,
    inc_angles,
    phase_names,
    color,
    height=None,
    horizontal=False,
):
    if horizontal:
        width = 15.0
        height = 6.0

        axis_height = 5.0 / height
        resid_heigt = 1.0 - axis_height
        title_height = resid_heigt

        axis_width = 5.0 / width
    else:
        if height == None:
            height = 19.0
        axis_height = 5.0 / height
        resid_height = 1.0 - 3.0 * axis_height
        title_height = resid_height / 3.0

    ## Full moment tensor:
    img1, buf1 = Get_bb_img(MT_Full, color)

    if horizontal:
        fig = plt.figure(figsize=(width, height), dpi=200)
        ax_1 = fig.add_axes([0.0, 0.0, axis_width, axis_height])
    else:
        fig = plt.figure(figsize=(5, height), dpi=200)
        ax_1 = fig.add_axes([0.0, 2 * (axis_height + title_height), 1.0, axis_height])

    if img1 is not None:
        ax_1.imshow(img1 / np.max(img1.flatten()))
    if horizontal:
        ax_X = fig.add_axes([0.0, 0.0, axis_width, axis_height], label="Circle_ray")
    else:
        ax_X = fig.add_axes(
            [0.0, 2 * (axis_height + title_height), 1.0, axis_height], label="Circle_ray"
        )
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)
    if azimuths is not None and inc_angles is not None:
        for a, i, phase in zip(azimuths, inc_angles, phase_names):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            p = Circle(
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=17)
    for a in [ax_1, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")
    #
    if horizontal:
        title_1 = fig.add_axes([0.0, axis_height, axis_width, title_height])
    else:
        title_1 = fig.add_axes([0.0, 3 * axis_height + 2 * title_height, 1.0, title_height])
    title_1.set_xticks([])
    title_1.set_yticks([])
    title_1.axis("off")
    title_1.text(
        0.5,
        0.2,
        "Full moment\n" r"$\epsilon=%.2f$" % Eps,
        ha="center",
        va="bottom",
        size="x-large",
        fontsize=25,
    )

    ########################

    ## DC moment tensor:
    img2, buf2 = Get_bb_img(MT_DC, color)
    if horizontal:
        ax_2 = fig.add_axes([axis_width, 0.0, axis_width, axis_height])
    else:
        ax_2 = fig.add_axes([0.0, axis_height + title_height, 1.0, axis_height])

    if img2 is not None:
        ax_2.imshow(img2 / np.max(img2.flatten()))
    if horizontal:
        ax_X = fig.add_axes([axis_width, 0.0, axis_width, axis_height], label="Circle_ray")
    else:
        ax_X = fig.add_axes(
            [0.0, axis_height + title_height, 1.0, axis_height], label="Circle_ray"
        )
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)
    if azimuths is not None and inc_angles is not None:
        for a, i, phase in zip(azimuths, inc_angles, phase_names):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            p = Circle(
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=17)
    for a in [ax_2, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")

    if horizontal:
        title_2 = fig.add_axes([axis_width, axis_height, axis_width, title_height])
    else:
        title_2 = fig.add_axes([0.0, 2 * axis_height + title_height, 1.0, title_height])
    title_2.set_xticks([])
    title_2.set_yticks([])
    title_2.axis("off")
    title_2.text(
        0.5,
        0.2,
        "Double-Couple \n M0: %.2e" % M0_DC,
        ha="center",
        va="bottom",
        size="x-large",
        fontsize=25,
    )
    # title_2.text(0.5, 0.2, "Direct", ha="center", va="bottom", size="x-large", fontsize=40)

    ### CLVD
    img3, buf3 = Get_bb_img(MT_CLVD, color)
    if horizontal:
        ax_3 = fig.add_axes([2 * (axis_width), 0.0, axis_width, axis_height])
    else:
        ax_3 = fig.add_axes([0.0, 0.0, 1.0, axis_height])

    if img3 is not None:
        ax_3.imshow(img3 / np.max(img3.flatten()))
    if horizontal:
        ax_X = fig.add_axes([2 * (axis_width), 0.0, axis_width, axis_height], label="Circle_ray")
    else:
        ax_X = fig.add_axes([0.0, 0.0, 1.0, axis_height], label="Circle_ray")
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)
    if azimuths is not None and inc_angles is not None:
        for a, i, phase in zip(azimuths, inc_angles, phase_names):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            p = Circle(
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=17)
    for a in [ax_3, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")

    if horizontal:
        title_3 = fig.add_axes([2 * axis_width, axis_height, axis_width, title_height])
    else:
        title_3 = fig.add_axes([0.0, axis_height, 1.0, title_height])
    title_3.set_xticks([])
    title_3.set_yticks([])
    title_3.axis("off")
    title_3.text(
        0.5,
        0.2,
        "CLVD \n M0: %.2e" % M0_CLVD,
        ha="center",
        va="bottom",
        size="x-large",
        fontsize=25,
    )

    return fig


def Plot_GS_BB(
    strikes, dips, rakes, azimuths, inc_angles, phase_names, color, height=None, horizontal=True
):
    if horizontal:
        width = 5.0
        height = 6.0

        axis_height = 5.0 / height
        resid_heigt = 1.0 - axis_height
        title_height = resid_heigt

        axis_width = 5.0 / width
    else:
        if height == None:
            height = 19.0
        axis_height = 5.0 / height
        resid_height = 1.0 - 3.0 * axis_height
        title_height = resid_height / 3.0

    fig_bb = plt.figure(figsize=(5, 5), dpi=200)
    ax_bb = fig_bb.add_axes([0.0, 0.0, 1.0, 1.0])

    ax_bb.set_xticks([])
    ax_bb.set_yticks([])
    ax_bb.axis("off")
    img = None
    buf = io.BytesIO()
    i = 0
    for strike, dip, rake in zip(strikes, dips, rakes):
        i += 1
        b = beach(
            fm=[strike, dip, rake],
            width=990,
            linewidth=0,
            facecolor=color,
            xy=(0, 0),
            axes=ax_bb,
            alpha=1,
            zorder=i,
        )
        ax_bb.add_collection(b)
        ax_bb.set_xlim((-1, 1))
        ax_bb.set_ylim((-1, 1))

        buf.seek(0)
        fig_bb.savefig(buf, format="png", dpi=200)
        buf.seek(0)
        if img is None:
            img = mpimg.imread(buf)
        else:
            img += mpimg.imread(buf)

    plt.close(fig_bb)

    if horizontal:
        fig = plt.figure(figsize=(width, height), dpi=200)
        ax_1 = fig.add_axes([0.0, 0.0, axis_width, axis_height])
    else:
        fig = plt.figure(figsize=(5, height), dpi=200)
        ax_1 = fig.add_axes([0.0, 2 * (axis_height + title_height), 1.0, axis_height])

    if img is not None:
        ax_1.imshow(img / np.max(img.flatten()))
    if horizontal:
        ax_X = fig.add_axes([0.0, 0.0, axis_width, axis_height], label="Circle_ray")
    else:
        ax_X = fig.add_axes(
            [0.0, 2 * (axis_height + title_height), 1.0, axis_height], label="Circle_ray"
        )
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)
    if azimuths is not None and inc_angles is not None:
        for a, i, phase in zip(azimuths, inc_angles, phase_names):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            p = Circle(
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=17)
    for a in [ax_1, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")
    #
    if horizontal:
        title_1 = fig.add_axes([0.0, axis_height, axis_width, title_height])
    else:
        title_1 = fig.add_axes([0.0, 3 * axis_height + 2 * title_height, 1.0, title_height])
    title_1.set_xticks([])
    title_1.set_yticks([])
    title_1.axis("off")
    title_1.text(0.5, 0.2, "Grid-Search", ha="center", va="bottom", size="x-large", fontsize=40)
    return fig

