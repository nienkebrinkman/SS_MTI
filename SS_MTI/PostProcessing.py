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
import glob
from os.path import join as pjoin
from typing import List as _List, Union as _Union
import instaseis
from obspy import UTCDateTime as utct
from obspy.imaging.beachball import aux_plane
from sklearn.cluster import KMeans
from pyrocko import moment_tensor as mtm
import matplotlib.patches as mpatches
import glob

pyproj_datadir = os.environ["PROJ_LIB"]

from mpl_toolkits.basemap import Basemap
import re

from SS_MTI import Read_H5 as _ReadH5
from SS_MTI import MTDecompose as _MTDecompose
from SS_MTI import Forward as _Forward
from SS_MTI import PreProcess as _PreProcess
from SS_MTI import GreensFunctions as _GreensFunctions
from SS_MTI import RadiationPattern as _RadiationPattern


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
                        event_vp, depth_syn[i], "r*", markersize=15, label="Synthetic Depth",
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
                        event_vs, depth_syn[i], "r*", markersize=15, label="Synthetic Depth",
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
                        event_dens, depth_syn[i], "r*", markersize=15, label="Synthetic Depth",
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
            handles=unique_list, prop={"size": 6}, loc="upper left", bbox_to_anchor=(0.0, 1.07),
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
        nrows=len(stream), ncols=len(phase_cuts), figsize=(18, 8), sharex="col", sharey="all",
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
                0, y, phase_cuts[j], verticalalignment="center", color="grey", fontsize=6,
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
        handles=unique_list, prop={"size": 6}, loc="upper left", bbox_to_anchor=(0.0, 1.4),
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
        projection="merc", llcrnrlat=-80, urcrnrlat=80, llcrnrlon=0, urcrnrlon=200, resolution="c",
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
    # EQlonB, EQlatB = m(lo_sB, la_sB) # 235b
    # EQlonC, EQlatC = m(lo_sC, la_sC)
    # EQlonD, EQlatD = m(lo_sC, la_sC)
    m.plot(EQlonA, EQlatA, "r*", markersize=20, zorder=10, label=name)
    # m.plot(EQlonB, EQlatB, 'g*', markersize=20, zorder=10, label = event_B.name)
    #     m.plot(EQlonC, EQlatC, 'b*', markersize=20, zorder=10, label=event_C.name)
    plt.legend(fontsize=20)
    plt.tight_layout()
    # plt.show()
    #     plt.savefig('Location_Event.pdf')
    return fig


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
            fm=MT, width=990, linewidth=0, facecolor=color, xy=(0, 0), axes=ax_bb_1, alpha=alpha,
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

    DC_scal = np.sqrt(1 - Eps / 0.5)
    CLVD_scal = np.sqrt(1 - (1 - Eps / 0.5))

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
            [0.0, 2 * (axis_height + title_height), 1.0, axis_height], label="Circle_ray",
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
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True,
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=40)
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
    # title_1.text(
    #     0.5,
    #     0.2,
    #     "Full moment\n" r"$\epsilon=%.2f$" % Eps,
    #     ha="center",
    #     va="bottom",
    #     size="x-large",
    #     fontsize=50,
    # )
    title_1.text(
        0.5, 0.2, "$\epsilon=%.2f$" % Eps, ha="center", va="bottom", size="x-large", fontsize=50,
    )

    ########################

    ## DC moment tensor:
    img2, buf2 = Get_bb_img(MT_DC, color)
    if horizontal:
        ax_2 = fig.add_axes(
            [
                axis_width + ((axis_width - (axis_width * DC_scal)) / 2),
                0.0 + ((axis_height - (axis_height * DC_scal)) / 2),
                axis_width * DC_scal,
                axis_height * DC_scal,
            ]
        )
    else:
        ax_2 = fig.add_axes([0.0, axis_height + title_height, 1.0, axis_height])

    if img2 is not None:
        ax_2.imshow(img2 / np.max(img2.flatten()))
    if horizontal:
        ax_X = fig.add_axes(
            [
                axis_width + ((axis_width - (axis_width * DC_scal)) / 2),
                0.0 + ((axis_height - (axis_height * DC_scal)) / 2),
                axis_width * DC_scal,
                axis_height * DC_scal,
            ],
            label="Circle_ray",
        )
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
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True,
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=40)
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
    # title_2.text(
    #     0.5,
    #     0.2,
    #     "Double-Couple \n M0: %.2e" % M0_DC,
    #     ha="center",
    #     va="bottom",
    #     size="x-large",
    #     fontsize=25,
    # )
    # title_2.text(0.5, 0.2, "Direct", ha="center", va="bottom", size="x-large", fontsize=40)

    ### CLVD
    img3, buf3 = Get_bb_img(MT_CLVD, color)
    if horizontal:
        ax_3 = fig.add_axes(
            [
                2 * (axis_width) + ((axis_width - (axis_width * CLVD_scal)) / 2),
                0.0 + ((axis_height - (axis_height * CLVD_scal)) / 2),
                axis_width * CLVD_scal,
                axis_height * CLVD_scal,
            ]
        )
    else:
        ax_3 = fig.add_axes([0.0, 0.0, 1.0, axis_height])

    if img3 is not None:
        ax_3.imshow(img3 / np.max(img3.flatten()))
    if horizontal:
        ax_X = fig.add_axes(
            [
                2 * (axis_width) + ((axis_width - (axis_width * CLVD_scal)) / 2),
                0.0 + ((axis_height - (axis_height * CLVD_scal)) / 2),
                axis_width * CLVD_scal,
                axis_height * CLVD_scal,
            ],
            label="Circle_ray",
        )
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
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True,
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=40)
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
    # title_3.text(
    #     0.5,
    #     0.2,
    #     "CLVD \n M0: %.2e" % M0_CLVD,
    #     ha="center",
    #     va="bottom",
    #     size="x-large",
    #     fontsize=25,
    # )

    return fig


def Plot_GS_BB(
    strikes, dips, rakes, azimuths, inc_angles, phase_names, color, height=None, horizontal=True,
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
            [0.0, 2 * (axis_height + title_height), 1.0, axis_height], label="Circle_ray",
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
                (x, y), 0.015, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True,
            )
            ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=40)
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
    # title_1.text(0.5, 0.2, "Grid-Search", ha="center", va="bottom", size="x-large", fontsize=40)
    return fig


""" Misfit analysis """


def plot_misfit_vs_depth(
    baz: float,
    save_paths: [str] = [],
    depths: [int] = [45],
    DOF: float = 700,
    event_name: str = "S0235b",
    misfit_name: str = "L2",
    veloc_model: str = "TAYAK_BKE",
    true_depth: float = None,
    Moho: float = 30,
    fmin: float = 1.0 / 10.0,
    fmax: float = 1.0 / 2.0,
    amount_of_phases: int = 5,
):
    labels = ["", ""]
    n_lowest = 1

    fig, ax = plt.subplots(
        nrows=3, ncols=1, sharex="all", figsize=(28, 17), gridspec_kw={"height_ratios": [4, 1, 1]},
    )
    # if event_name == "S0183a":
    #     fig, ax = plt.subplots(nrows=1, ncols=1, sharex="all", figsize=(28, 5.33),)
    # else:
    #     fig, ax = plt.subplots(
    #         nrows=3,
    #         ncols=1,
    #         sharex="all",
    #         figsize=(28, 8),
    #         gridspec_kw={"height_ratios": [4, 1, 1]},
    #     )
    # from matplotlib import gridspec
    # gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    BB = []
    Line_x = []
    Line1_ymin = []
    Line1_ymax = []
    Line2_ymin = []
    Line2_ymax = []
    for i, save_path in enumerate(save_paths):
        L2_GS = np.array([])
        L2_Direct = np.array([])
        Eps = np.array([])
        cond_nrs = np.array([])

        for idepth, depth in enumerate(depths):
            print(i, depth)
            GS_File = glob.glob(
                pjoin(
                    save_path,
                    f"GS_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_model}.hdf5",
                )
            )[0]
            if event_name == "S0183a":
                pass
            else:
                Direct_File = glob.glob(
                    pjoin(
                        save_path,
                        f"Direct_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{veloc_model}.hdf5",
                    )
                )[0]

            ## ================ READ GS =============================
            (depth_GS, sdr, M0_GS, misfit_L2_GS,) = _ReadH5.Read_GS_h5(
                Filename=GS_File, amount_of_phases=amount_of_phases
            )
            Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
            # Total_L2_norm_GS = np.sum(misfit_L2_norm_GS, axis=1)
            # GOF = ( (Total_L2_GS - DOF ) * 100 ) / DOF
            # GOF_GS = (Total_L2_norm_GS / DOF) * 100
            GOF_GS = Total_L2_GS / DOF

            lowest_indices = Total_L2_GS.argsort()[0:n_lowest]
            # lowest_indices = GOF_GS.argsort()[0:n_lowest]

            sdr = sdr[lowest_indices, :]
            print("strike", sdr[0][0], "dip", sdr[0][1], "rake", sdr[0][2])
            depth_GS = depth_GS[lowest_indices]
            M0_GS = M0_GS[lowest_indices]
            # shifts["P"] = shifts["P"][lowest_indices]
            # shifts["S"] = shifts["S"][lowest_indices]

            # L2_GS = np.append(L2_GS, Total_L2_GS[lowest_indices][0])
            L2_GS = np.append(L2_GS, GOF_GS[lowest_indices][0])
            # if depth == 8:
            #     lowest_indices[0] = lowest_indices[2]
            #     L2_GS = np.append(L2_GS, GOF_GS[lowest_indices][0])
            # else:
            #     L2_GS = np.append(L2_GS, GOF_GS[lowest_indices][0])
            if event_name == "S0183a":
                pass
            else:
                ## ============ Read Direct ========================
                (
                    depth_Direct,
                    FULL_MT,
                    DC_MT,
                    CLVD_MT,
                    misfit_L2_Direct,
                    Epsilon,
                    M0,
                    M0_DC,
                    M0_CLVD,
                    angles,
                    cond_nr,
                ) = _ReadH5.Read_Direct_Inversion(Direct_File, amount_of_phases=amount_of_phases)

                Total_L2_Direct = np.sum(misfit_L2_Direct)
                # Total_L2_norm_Direct = np.sum(misfit_L2_norm_Direct, axis=1)
                # GOF_Direct = (Total_L2_norm_Direct / DOF) * 100
                GOF_Direct = Total_L2_Direct / DOF
                # L2_Direct = np.append(L2_Direct, Total_L2_Direct[0])
                L2_Direct = np.append(L2_Direct, GOF_Direct)

            # ============== CREATE BEACHBALL PATCHES ===============

            # y1 = Total_L2_GS[lowest_indices][0]
            # y2 = Total_L2_Direct[0]

            y2 = GOF_GS[lowest_indices][0]
            if event_name == "S0183a":
                ax_current = ax
                pass
            else:
                y1 = GOF_Direct
                ax_current = ax[0]
                y_dist = np.log(np.abs(y1 - y2))
                if y_dist < 100:
                    adding_value = 2e-1
                    Line_x.append(depth)
                    if y1 > y2:
                        y1 = y1 + adding_value
                        y2 = y2 - adding_value
                        Line1_ymin.append(GOF_Direct)
                        Line1_ymax.append(y1)
                        Line2_ymin.append(y2)
                        Line2_ymax.append(GOF_GS[lowest_indices][0])
                    else:
                        diff = y2 - y1
                        y1 = y1 + adding_value + diff
                        y2 = y2 - adding_value - diff
                        Line1_ymax.append(GOF_Direct)
                        Line1_ymin.append(y1)
                        Line2_ymax.append(y2)
                        Line2_ymin.append(GOF_GS[lowest_indices][0])

            BB.append(
                beach(
                    [sdr[0][0], sdr[0][1], sdr[0][2]],
                    xy=(depth_GS[0], y2),
                    width=40,
                    linewidth=1,
                    axes=ax_current,
                )
            )
            if event_name == "S0183a":
                pass
            else:
                BB.append(
                    beach(
                        DC_MT / M0_DC,
                        xy=(depth, y1),
                        width=40,
                        facecolor="r",
                        linewidth=1,
                        axes=ax_current,
                    )
                )

                Eps = np.append(Eps, Epsilon)
                cond_nrs = np.append(cond_nrs, cond_nr)

        ax_current.plot(depths, L2_GS, "-bo", label="Grid-Search %s" % labels[i], lw=i + 1)
        if event_name == "S0183a":
            pass
        else:
            ax_current.plot(depths, L2_Direct, "-ro", label="Direct %s" % labels[i], lw=i + 1)
        if i == 0:
            ax_current.axvline(x=Moho, c="grey", ls="dashed", label="Moho", lw=3)
            # true_depth = 45.
            if true_depth is not None:
                ax_current.axvline(x=true_depth, c="green", ls="dotted", label="True Depth", lw=2)

        if event_name == "S0183a":
            pass
        else:
            ax[1].plot(depths, Eps, "--ko", label="Epsilon %s" % labels[i], lw=0.5)
            if i == 0:
                ax[1].axvline(x=Moho, c="grey", ls="dashed", lw=3)
                if true_depth is not None:
                    ax[1].axvline(x=true_depth, c="green", ls="dotted", label="True Depth")
            # ax[2].semilogy(depths, cond_nrs, "--ko", label="Condition number %s" % labels[i], lw=0.5)
            ax[2].plot(depths, cond_nrs, "--ko", label="Condition number %s" % labels[i], lw=0.5)
            ax[2].ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
            # if event_name == "S0235b":
            #     ax[2].set_yticks([80, 100, 200, 400])
            # elif event_name == "S0173a":
            #     ax[2].set_yticks([700, 1000, 1300])

            if i == 0:
                ax[2].axvline(x=Moho, c="grey", ls="dashed", lw=3)
                if true_depth is not None:
                    ax[2].axvline(x=true_depth, c="green", ls="dotted", label="True Depth")

        for iline in range(len(Line_x)):
            ax_current.plot(
                [Line_x[iline], Line_x[iline]],
                [Line1_ymin[iline], Line1_ymax[iline]],
                c="r",
                ls="dashed",
                alpha=0.5,
                lw=0.5,
            )
            ax_current.plot(
                [Line_x[iline], Line_x[iline]],
                [Line2_ymin[iline], Line2_ymax[iline]],
                c="b",
                ls="dashed",
                alpha=0.5,
                lw=0.5,
            )
    for bb in BB:
        ax_current.add_collection(bb)

    ax_current.legend(prop={"size": 45}, loc="upper center", ncol=len(save_paths) + 1)
    ax_current.set_ylabel(r"$\chi^2$", fontsize=45)
    ax_current.tick_params(axis="both", which="major", labelsize=23)
    ax_current.tick_params(axis="both", which="minor", labelsize=23)
    ax_current.grid(True)
    # ax_current.set_ylim(-0.5, 10)

    if event_name == "S0183a":
        pass
    else:
        extraticks = [0.1, 0.2, 0.3, 0.4]
        ax[1].set_yticks(list(ax[1].get_yticks()) + extraticks)
        # ax[1].legend(prop={"size": 15}, loc="upper right")
        # ax[1].set_xlabel("Depth (km)", fontsize=20)
        ax[1].set_ylabel(r"$\epsilon$", fontsize=45)
        ax[1].tick_params(axis="both", which="major", labelsize=23)
        ax[1].tick_params(axis="both", which="minor", labelsize=23)
        ax[1].set_ylim(-0.05, 0.5)
        ax[1].grid(True)

        ax[2].set_xlabel("Depth (km)", fontsize=45)
        ax[2].set_ylabel(r"$\kappa$", fontsize=45)
        ax[2].tick_params(axis="both", which="major", labelsize=23)
        ax[2].tick_params(axis="both", which="minor", labelsize=23)
        # if event_name == "S0173a":
        #     ax[2].set_ylim(0.0, 2000.0)
        ax[2].grid(True)
    return fig


def plot_phases_vs_depth(
    h5_file_folder: str,
    method: str,
    misfit_name: str,
    fwd: _Forward._AbstractForward,
    event: obspy.core.event.Event,
    rec: instaseis.Receiver,
    phases: [str],
    components: [str],
    t_pre: [float],
    t_post: [float],
    depths: [float],
    phase_corrs: [float] = None,
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = None,
    tstars: _Union[_List[float], _List[str]] = None,
    color_plot: str = None,
    pref_depth_start=[None],
    pref_depth_end=[None],
):
    assert method == "GS" or method == "Direct", "method needs to be either GS or Direct"

    if tstars is None:
        tstars = [None] * len(phases)

    if phase_corrs is None:
        phase_corrs = [0] * len(phases)

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    # TODO: IMPLEMENT LQT COORDINATE SYSTEM
    LQT_value = False
    baz = None
    inc = None

    SHIFT = 2.5
    """ Process the observed data """
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
    t_post_new = [t + SHIFT for t in t_post]
    st_obs, sigmas = _PreProcess.prepare_event_data(
        event=event,
        phases=phases,
        components=components,
        slice=True,
        tts=obs_tt,
        t_pre=t_pre,
        t_post=t_post_new,
        filter=filter_par,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        noise_level=False,
    )

    widths = [5] * len(phases) + [1]
    fig, ax = plt.subplots(
        nrows=1,
        ncols=len(phases) + 1,
        figsize=(4 * len(phases), 4 * (len(phases) + 1)),
        sharex="col",
        sharey="col",
        gridspec_kw={"width_ratios": widths},
    )
    Yticks = np.arange(len(depths) + len(pref_depth_start)) * 1.8
    obs = 0
    pref_depth_end_sorted = pref_depth_end.sort()
    pref_depth_start.sort()
    BB = []
    extra_arrs = [[] for _ in range(len(depths))]
    syn_tts = [[] for _ in range(len(depths))]
    for idepth, depth in enumerate(depths):
        print(depth)
        if method == "GS":
            if event.name == "S0173a":
                MT_depth = 38
            elif event.name == "S0235b":
                MT_depth = 32  # 14
            h5_file_path = pjoin(
                h5_file_folder,
                f"GS_{event.name}_{MT_depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )
            depth_GS, sdr, M0_GS, misfit_L2_GS = _ReadH5.Read_GS_h5(
                Filename=h5_file_path, amount_of_phases=5
            )
            Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
            n_lowest = 1
            lowest_indices = Total_L2_GS.argsort()[0:n_lowest]
            MT = sdr[lowest_indices, :][0]
            # MT = [340.0, 90.0, 105.0]
            print("strike", MT[0], "dip", MT[1], "rake", MT[2])
            depth_GS = depth_GS[lowest_indices]
            M0 = M0_GS[lowest_indices][0]
        else:
            h5_file_path = pjoin(
                h5_file_folder,
                f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )
            (
                depth_Direct,
                MT,
                DC_MT,
                CLVD_MT,
                misfit_L2_Direct,
                Epsilon,
                M0,
                M0_DC,
                M0_CLVD,
                angles,
                cond_nr,
            ) = _ReadH5.Read_Direct_Inversion(h5_file_path, amount_of_phases=5)
            Total_L2_Direct = np.sum(misfit_L2_Direct)

            MT = DC_MT

        """ Generate Green's functions per depth """

        Underside_refl_src = []
        Conversion_src = []
        Conversion_rec = []
        Reflection_phases = []
        model = fwd.taup_veloc
        for i, values in enumerate(model.model.s_mod.critical_depths):
            if values[0] < 1.0 or values[0] > 100.0:
                continue
            interface = str(int(values[0]))
            if i > 1:
                Reflection_phases.append(
                    "P"
                    + interface
                    + "s"
                    + str(int(model.model.s_mod.critical_depths[i - 1][0]))
                    + "p"
                )
            # if values[0] > 50.0 and "TAYAK" in model_name:
            #     continue
            for down_phase in ["p^", "s^"]:
                for up_phase in ["P", "S"]:

                    Underside_refl_src.append(down_phase + interface + up_phase)

            Conversion_src.append("S" + interface + "P")
            Conversion_src.append("P" + interface + "S")
            Conversion_rec.append("P" + interface + "p")
            Conversion_rec.append("P" + interface + "s")
        # extra_phases = Conversion_src
        if event.name == "S0173a":
            add = ["p^24P", "p^10P", "s^24S", "s^10S", "P10s", "P24s"]
        elif event.name == "S0235b":
            add = ["p^24P", "p^10P", "s^24S", "s^10S", "P10s", "P24s"]
        extra_phases = ["pP", "sS", "sP", "pS", "SS", "PP", "SSS", "PPP",] + add
        if not os.path.exists(os.path.join(h5_file_folder, "ray_paths")):
            os.mkdir(os.path.join(h5_file_folder, "ray_paths"))
        for j, extraphase in enumerate(extra_phases):
            arrivals = fwd.taup_veloc.get_ray_paths(
                source_depth_in_km=depth,
                distance_in_degree=event.distance,
                phase_list=[extraphase],
            )
            if arrivals:
                ax_ray = arrivals.plot_rays(plot_type="cartesian", show=False, legend=True)
                plt.savefig(os.path.join(h5_file_folder, "ray_paths", f"d_{depth}_{extraphase}"))
                plt.close()
            extra_arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
            if extra_arr:
                extra_arrs[idepth].append(extra_arr)
            else:
                extra_arrs[idepth].append(extra_arr)

        for i, phase in enumerate(phases):
            syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)
            syn_tts[idepth].append(syn_tt)
            syn_GF = fwd.get_greens_functions(
                comp=components[i],
                depth=depth,
                distance=event.distance,
                lat_src=event.latitude,
                lon_src=event.longitude,
                rec=rec,
                tstar=tstars[i],
                LQT=LQT_value,
                inc=inc,
                baz=baz,
                M0=1.0,
                filter=filter_par,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
            )
            tr_syn = fwd.generate_synthetic_data(
                st_GF=syn_GF,
                focal_mech=MT,
                M0=M0,
                slice=True,
                tt=syn_tt,
                t_pre=t_pre[i],
                t_post=t_post[i] + SHIFT,
            )

            color = "k"
            lw = 1
            if (
                all(x is None for x in pref_depth_start) is not None
                and all(x is None for x in pref_depth_end) is not None
            ):
                for d_range in range(len(pref_depth_end)):
                    if depth <= pref_depth_end[d_range] and depth >= pref_depth_start[d_range]:
                        color = "purple"
                        lw = 2

            ytick = Yticks[idepth]
            if (depth == depths[[np.abs(depths - x).argmin() + 1 for x in pref_depth_end]]).any():
                obs = (
                    1
                    + np.where(
                        depth == depths[[np.abs(depths - x).argmin() + 1 for x in pref_depth_end]]
                    )[0][0]
                )
                print(obs)
                ytick = Yticks[idepth + obs - 1]
                fig, ax[i] = Plot_phase_vs_depth_copy(
                    tr=st_obs[i],
                    depth=depth,
                    total_depths=len(depths) + len(pref_depth_start),
                    Ytick=ytick,
                    phase=phase,
                    t_pre=t_pre[i],
                    t_post=t_post[i],
                    fig=fig,
                    ax=ax[i],
                    extra_phases=None,
                    extra_arrs=None,
                    color="b",
                    linewidth=2,
                    SHIFT=SHIFT,
                )
                ytick = Yticks[idepth + obs]
            elif (depth > depths[[np.abs(depths - x).argmin() + 1 for x in pref_depth_end]]).any():
                ytick = Yticks[idepth + obs]

            fig, ax[i] = Plot_phase_vs_depth_copy(
                tr=tr_syn,
                depth=depth,
                total_depths=len(depths) + len(pref_depth_start),
                Ytick=ytick,
                phase=phase,
                t_pre=t_pre[i],
                t_post=t_post[i],
                fig=fig,
                ax=ax[i],
                extra_phases=extra_phases,
                extra_arrs=extra_arrs[idepth],
                color=color,
                linewidth=lw,
                SHIFT=SHIFT,
            )
            BB.append(
                beach(
                    MT,
                    xy=(0, ytick),
                    width=20,
                    linewidth=1,
                    alpha=0.5,
                    facecolor="r",
                    axes=ax[-1],
                )
            )

    delta = Yticks[1] - Yticks[0]

    idxs = [np.abs(depths - x).argmin() + 1 for x in pref_depth_end]

    depths_ins = depths
    for idx_i, idx in enumerate(idxs):
        depths_ins = np.insert(depths_ins, idx + idx_i, 70)
    for idx_i, idx in enumerate(idxs):
        for i, phase in enumerate(phases):
            if idx_i > 0:
                fill_zero = idxs[idx_i]
            else:
                fill_zero = 0
            if idx_i == len(idxs) - 1:
                fill_one = -1
                fill_zero = idxs[idx_i] + len(idxs)
            else:
                fill_one = idxs[idx_i + 1]
            if phase == "P":
                ymax = ax[i].get_ylim()[1]
                if event.name == "S0173a":
                    if idx_i == 0:
                        ax[i].fill(
                            [17 - SHIFT, 40 - SHIFT, 40 - SHIFT, 17 - SHIFT],
                            [
                                Yticks[fill_zero] - 0.4 * delta,
                                Yticks[fill_zero] - 0.4 * delta,
                                Yticks[idx - (idx_i + 1)] + 0.4 * delta,
                                Yticks[idx - (idx_i + 1)] + 0.4 * delta,
                            ],
                            facecolor="green",
                            alpha=0.2,
                        )
                    ax[i].fill(
                        [17 - SHIFT, 40 - SHIFT, 40 - SHIFT, 17 - SHIFT],
                        [
                            Yticks[idx + (idx_i + 1)] - 0.4 * delta,
                            Yticks[idx + (idx_i + 1)] - 0.4 * delta,
                            Yticks[fill_one] + 0.4 * delta,
                            Yticks[fill_one] + 0.4 * delta,
                        ],
                        facecolor="green",
                        alpha=0.2,
                    )
                    ax[i].text(
                        17,
                        ymax * 0.8,
                        "Glitch",
                        verticalalignment="center",
                        color="green",
                        fontsize=8,
                    )

            # fig, ax[i] = Plot_phase_vs_depth_copy(
            #     tr=st_obs[i],
            #     depth=depth,
            #     total_depths=len(depths) + 1,
            #     Ytick=Yticks[idepth + 1],
            #     phase=phase,
            #     t_pre=t_pre[i],
            #     t_post=t_post[i],
            #     fig=fig,
            #     ax=ax[i],
            #     extra_phases=None,
            #     extra_arrs=None,
            #     color="b",
            # )
            # ax[i].axvline(0.0, c="dimgrey")
            if idx_i == 0:
                ax[i].plot(
                    [0, 0],
                    [Yticks[fill_zero] - 0.4 * delta, Yticks[idx - 1] + 0.4 * delta],
                    "dimgrey",
                    lw=1,
                )
                ax[i].text(
                    0.1,
                    Yticks[idx - 1] + 0.4 * delta,
                    phase,
                    verticalalignment="center",
                    color="dimgrey",
                    fontsize=15,
                    weight="bold",
                )
            ax[i].plot(
                [0, 0],
                [Yticks[idx + (idx_i + 1)] - 0.4 * delta, Yticks[fill_one] + 0.4 * delta],
                "dimgrey",
                lw=1,
            )
            ax[i].text(
                0.1,
                Yticks[-1] + 0.4 * delta,
                phase,
                verticalalignment="center",
                color="dimgrey",
                fontsize=15,
                weight="bold",
            )

            for k in range(len(extra_phases)):
                """ Below observed plot: """
                if idx_i == 0:
                    syn_tt_depth = np.asarray(
                        [arr[i] for arr in syn_tts[fill_zero:idx]], dtype=np.float
                    )
                    # syn_tt_depth = np.asarray([arr[i] for arr in syn_tts], dtype=np.float)
                    x = np.asarray([arr[k] for arr in extra_arrs[fill_zero:idx]], dtype=np.float)
                    x = x - syn_tt_depth
                    if np.any(x > t_post[i]):
                        x[x > t_post[i]] = None
                    if np.any(x < -t_pre[i]):
                        x[x < -t_pre[i]] = None
                    # y = np.asarray(Yticks[:-1])
                    y = np.asarray(Yticks[fill_zero:idx])
                    y = y[~np.isnan(x)]
                    x = x[~np.isnan(x)]
                    if x.size == 0:
                        continue
                    rotn = np.degrees(np.arctan(y[-1:] - y[-2:-1], x[-1:] - x[-2:-1]))
                    if rotn.size == 0:
                        trans_angle = 90
                        ax[i].plot(
                            [x[0], x[0]], [y[0] - 0.4 * delta, y[0] + 0.4 * delta], "dimgrey", lw=1
                        )
                    else:
                        l2 = np.array((x[-1], y[-1]))
                        rotation = rotn[-1]
                        trans_angle = plt.gca().transData.transform_angles(
                            np.array((rotation,)), l2.reshape((1, 2))
                        )[0]
                        ax[i].plot(x, y, "-", c="dimgrey", lw=1)
                    ax[i].text(
                        x[-1],
                        y[-1],
                        extra_phases[k],
                        verticalalignment="center",
                        color="dimgrey",
                        fontsize=15,
                        rotation=trans_angle,
                        weight="bold",
                    )

                """ Above observed plot"""
                if idx_i == len(idxs) - 1:
                    syn_tt_depth = np.asarray([arr[i] for arr in syn_tts[idx:]], dtype=np.float)
                    # syn_tt_depth = np.asarray([arr[i] for arr in syn_tts], dtype=np.float)
                    x = np.asarray([arr[k] for arr in extra_arrs[idx:]], dtype=np.float)
                else:
                    syn_tt_depth = np.asarray(
                        [arr[i] for arr in syn_tts[idx:fill_one]], dtype=np.float
                    )
                    # syn_tt_depth = np.asarray([arr[i] for arr in syn_tts], dtype=np.float)
                    x = np.asarray([arr[k] for arr in extra_arrs[idx:fill_one]], dtype=np.float)
                x = x - syn_tt_depth
                if np.any(x > t_post[i]):
                    x[x > t_post[i]] = None
                if np.any(x < -t_pre[i]):
                    x[x < -t_pre[i]] = None
                print(idx_i)
                if idx_i == len(idxs) - 1:
                    y = np.asarray(Yticks[idx + (idx_i + 1) :])
                else:
                    y = np.asarray(Yticks[idx + (idx_i + 1) : fill_one + (idx_i + 1)])
                y = y[~np.isnan(x)]
                x = x[~np.isnan(x)]
                if x.size == 0:
                    continue
                rotn = np.degrees(np.arctan(y[-1:] - y[-2:-1], x[-1:] - x[-2:-1]))
                if rotn.size == 0:
                    trans_angle = 90
                    ax[i].plot(
                        [x[0], x[0]], [y[0] - 0.4 * delta, y[0] + 0.4 * delta], "dimgrey", lw=1
                    )
                else:
                    l2 = np.array((x[-1], y[-1]))
                    rotation = rotn[-1]
                    trans_angle = plt.gca().transData.transform_angles(
                        np.array((rotation,)), l2.reshape((1, 2))
                    )[0]
                    ax[i].plot(x, y, "-", c="dimgrey", lw=1)

                ax[i].text(
                    x[-1],
                    y[-1],
                    extra_phases[k],
                    verticalalignment="center",
                    color="dimgrey",
                    fontsize=15,
                    rotation=trans_angle,
                    weight="bold",
                )

                # depths_ins = np.insert(depths, idx + idx_i, 70)
                ax[i].yaxis.set_ticks(Yticks)
                ax[i].set_yticklabels(depths_ins)
                yticks = ax[i].yaxis.get_major_ticks()
                yticks[idx + idx_i].set_visible(False)
                ax[i].set_ylim(Yticks[0] - delta, Yticks[-1] + delta)
                if not i == 0:
                    ax[i].get_yaxis().set_visible(False)

    # # for i in range(len(phases)):

    ax[-1].yaxis.set_ticks(Yticks)
    ax[-1].set_yticklabels(depths_ins)
    # yticks = ax[-1].yaxis.get_major_ticks()
    # yticks[idx].set_visible(False)
    ax[-1].set_ylim(Yticks[0] - delta, Yticks[-1] + delta)
    ax[-1].axis("off")
    ax[-1].set_xlim(-0.15, 0.15)
    for bb in BB:
        ax[-1].add_collection(bb)
    fig.text(
        0.5,
        0.95,
        f"Event {event.name}",
        ha="center",
        va="bottom",
        size="x-large",
        color="blue",
        fontsize=25,
    )
    fig.text(0.04, 0.5, "Source Depth (km)", va="center", rotation="vertical", fontsize=25)
    fig.text(0.5, 0.04, "Time after arrival (s)", va="center", fontsize=25)
    return fig


def Plot_phase_vs_depth_copy(
    tr: obspy.Trace,
    depth: float,
    total_depths: int,
    Ytick: float,
    phase: str,
    t_pre: float = 10.0,
    t_post: float = 50.0,
    fig: plt.figure = None,
    ax: plt.axes = None,
    extra_phases: [str] = None,
    extra_arrs: [float] = None,
    color: str = None,
    linewidth: float = 1,
    SHIFT: float = 0.0,
):
    if color is None:
        color = "k"
    st = tr.copy()
    # norm = st.slice(
    #     starttime=st.stats.starttime, endtime=st.stats.starttime + 7.0
    # ).max()  # TODO: this value needs to come from the function
    # st.data = st.data / norm
    st.normalize()
    global_max = st.data.max()
    global_min = st.data.min()
    y = global_max * 0.9 + Ytick
    ymin = global_min + Ytick
    ymax = global_max + Ytick

    ax.plot(st.times() - t_pre - SHIFT, st.data + Ytick, color, lw=linewidth)

    # ax.plot(
    #     [0, 0], [ymin, ymax], "grey",
    # )
    # ax.text(0, y, phase, verticalalignment="center", color="grey", fontsize=6)
    if extra_phases is not None:
        for k in range(len(extra_phases)):
            if extra_arrs[k] is None or extra_arrs[k] > t_post or extra_arrs[k] < -t_pre:
                continue
            phase_t = extra_arrs[k]

            y = global_max * 0.9 + Ytick
            c = "grey"

            ax.plot(
                [phase_t, phase_t], [ymin, ymax], c,
            )
            ax.text(
                phase_t + 0.1,
                y,
                extra_phases[k],
                verticalalignment="center",
                color=c,
                fontsize=6,
                rotation=90,
            )
    ax.set_xlim(-t_pre, t_post)
    ax.set_ylim(ymin, ymax)
    ax.set_title(f"{phase}-Phase channel:{st.stats.channel[-1]}")

    # fig.legend(handles=unique_list, prop={"size": 6}, loc="upper left")

    return fig, ax


def post_waveform_plotting(
    h5_file_folder: str,
    method: str,
    misfit_name: str,
    misfit_weight_len: float,
    fwd: _Forward._AbstractForward,
    event: obspy.core.event.Event,
    rec: instaseis.Receiver,
    phases: [str],
    components: [str],
    t_pre: [float],
    t_post: [float],
    depths: [float],
    phase_corrs: [float] = None,
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = None,
    tstars: _Union[_List[float], _List[str]] = None,
    plot_extra_phases: [str] = None,
    Ylims: [float] = None,
    Return_Fig: bool = False,
):
    if tstars is None:
        tstars = [None] * len(phases)

    if phase_corrs is None:
        phase_corrs = [0] * len(phases)

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    # TODO: IMPLEMENT LQT COORDINATE SYSTEM
    LQT_value = False
    baz = None
    inc = None

    """ PRE-PROCESS THE OBSERVED Travel_times """
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])

    for idepth, depth in enumerate(depths):
        print(depth)
        if method == "GS":
            color_plot = "b"
            h5_file_path = pjoin(
                h5_file_folder,
                f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )
            depth_GS, sdr, M0_GS, misfit_L2_GS = _ReadH5.Read_GS_h5(
                Filename=h5_file_path, amount_of_phases=len(phases)
            )
            Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
            lowest_ind = Total_L2_GS.argsort()
            Total_L2_GS.sort()
            misfit_low = Total_L2_GS[:] - Total_L2_GS[0]
            uncert = 0.05 * Total_L2_GS[0]
            inds = np.where(misfit_low < uncert)
            lowest_indices = lowest_ind[inds][0:1]
            print(np.sum(misfit_L2_GS, axis=1)[lowest_indices])
            # n_lowest = int(len(Total_L2_GS) * 0.05)
            # lowest_indices = Total_L2_GS.argsort()[0:n_lowest:50]
            # n_lowest = 10
            # lowest_indices = Total_L2_GS.argsort()[0:n_lowest]
            MT = sdr[lowest_indices, :]
            depth_GS = depth_GS[lowest_indices]
            M0 = M0_GS[lowest_indices]

            """ Calculate take-off angles"""
            takeoff_angles = ["P", "S", "pP"]
            angles = []
            for phase in takeoff_angles:
                angles.append(
                    fwd.get_phase_tt(
                        phase=phase, depth=depth, distance=event.distance, takeoffs=True
                    )
                )
            """ Beachball plot """
            fig = Plot_GS_BB(
                MT[:, 0],
                MT[:, 1],
                MT[:, 2],
                azimuths=[event.az, event.az, event.az],
                inc_angles=angles,
                phase_names=takeoff_angles,
                color=color_plot,
            )
            plt.savefig(
                pjoin(
                    h5_file_folder,
                    f"GS_BBB__{event.name}_{depth}_{misfit_name}_{fwd.veloc_name}_Post.pdf",
                ),
                dpi=300,
            )
            plt.close()
        else:
            h5_file_path = pjoin(
                h5_file_folder,
                f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )

            (
                depth_Direct,
                MT,
                DC_MT,
                CLVD_MT,
                misfit_L2_Direct,
                Epsilon,
                M0,
                M0_DC,
                M0_CLVD,
                angles,
                cond_nr,
            ) = _ReadH5.Read_Direct_Inversion(h5_file_path, amount_of_phases=5)
            Total_L2_Direct = np.sum(misfit_L2_Direct)

            color_plot = "r"

            MT_FULL_ = MT  # np.array([MT[0], MT[2], MT[1], MT[4], MT[3], MT[5],])
            DC_MT_ = (
                DC_MT  # np.array([DC_MT[0], DC_MT[2], DC_MT[1], DC_MT[4], DC_MT[3], DC_MT[5],])
            )
            CLVD_MT_ = CLVD_MT  # np.array(
            #     [CLVD_MT[0], CLVD_MT[2], CLVD_MT[1], CLVD_MT[4], CLVD_MT[3], CLVD_MT[5],]
            # )

            fig = Plot_Direct_BB(
                MT_Full=MT_FULL_ / M0,
                Eps=Epsilon,
                MT_DC=DC_MT_ / M0_DC,
                M0_DC=M0_DC,
                MT_CLVD=CLVD_MT_ / M0_CLVD,
                M0_CLVD=M0_CLVD,
                azimuths=[event.az, event.az, event.az],
                inc_angles=angles,
                phase_names=["P", "S", "pP"],
                color=color_plot,
                height=19.0,
                horizontal=True,
            )

            plt.savefig(
                pjoin(
                    h5_file_folder,
                    f"Direct_BB_{event.name}_{depth}_{misfit_name}_{fwd.veloc_name}_Post.pdf",
                ),
                dpi=300,
            )
            plt.close()

            MT = np.expand_dims(DC_MT, axis=0)
            M0 = np.expand_dims(M0_DC, axis=0)

        # """ Generate Green's functions per depth """
        # syn_tts = []
        # syn_GFs = []
        # for i, phase in enumerate(phases):
        #     syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)

        #     syn_GF = fwd.get_greens_functions(
        #         comp=components[i],
        #         depth=depth,
        #         distance=event.distance,
        #         lat_src=event.latitude,
        #         lon_src=event.longitude,
        #         rec=rec,
        #         tstar=tstars[i],
        #         LQT=LQT_value,
        #         inc=inc,
        #         baz=baz,
        #         M0=1.0,
        #         filter=filter_par,
        #         fmin=fmin,
        #         fmax=fmax,
        #         zerophase=zerophase,
        #     )
        #     syn_GFs.append(syn_GF)
        #     syn_tts.append(syn_tt)

        # if plot_extra_phases is not None:
        #     extra_arrs = []
        #     for j, extraphase in enumerate(plot_extra_phases):
        #         arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
        #         extra_arrs.append(arr)
        # else:
        #     extra_arrs = None

        # fig = waveform_plot(
        #     syn_GFs=syn_GFs,
        #     syn_tts=syn_tts,
        #     obs_tts=obs_tt,
        #     fwd=fwd,
        #     misfit_weight_len=misfit_weight_len,
        #     event=event,
        #     phases=phases,
        #     components=components,
        #     t_pre=t_pre,
        #     t_post=t_post,
        #     MTs=MT,
        #     M0s=M0,
        #     fmin=fmin,
        #     fmax=fmax,
        #     zerophase=zerophase,
        #     plot_extra_phases=plot_extra_phases,
        #     extra_arrs=extra_arrs,
        #     color_plot=color_plot,
        #     Ylims=Ylims,
        # )
        # if Return_Fig:
        #     return fig
        # else:
        #     plt.savefig(
        #         pjoin(
        #             h5_file_folder,
        #             f"{method}_waveforms_{event.name}_{depth}_{misfit_name}_{fwd.veloc_name}_Post.pdf",
        #         ),
        #         dpi=300,
        #     )
        #     plt.close()


def post_waveform_plotting_COMBINED(
    h5_file_folder: str,
    misfit_name: str,
    misfit_weight_len: float,
    fwd: _Forward._AbstractForward,
    event: obspy.core.event.Event,
    rec: instaseis.Receiver,
    phases: [str],
    components: [str],
    t_pre: [float],
    t_post: [float],
    depths: [float],
    phase_corrs: [float] = None,
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = None,
    tstars: _Union[_List[float], _List[str]] = None,
    plot_extra_phases: [str] = None,
    Ylims: [float] = None,
    Return_Fig: bool = False,
):
    if tstars is None:
        tstars = [None] * len(phases)

    if phase_corrs is None:
        phase_corrs = [0] * len(phases)

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    # TODO: IMPLEMENT LQT COORDINATE SYSTEM
    LQT_value = False
    baz = None
    inc = None

    """ PRE-PROCESS THE OBSERVED Travel_times """
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])

    for idepth, depth in enumerate(depths):
        print(depth)

        color_plot = "b"
        h5_file_path = pjoin(
            h5_file_folder,
            f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
        )
        depth_GS, sdr, M0_GS, misfit_L2_GS = _ReadH5.Read_GS_h5(
            Filename=h5_file_path, amount_of_phases=5
        )
        Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
        lowest_ind = Total_L2_GS.argsort()
        Total_L2_GS.sort()
        misfit_low = Total_L2_GS[:] - Total_L2_GS[0]
        uncert = 0.05 * Total_L2_GS[0]
        inds = np.where(misfit_low < uncert)
        lowest_indices = lowest_ind[inds][:50]

        # n_lowest = int(len(Total_L2_GS) * 0.05)
        # lowest_indices = Total_L2_GS.argsort()[0:n_lowest:50]
        # n_lowest = 10
        # lowest_indices = Total_L2_GS.argsort()[0:n_lowest]
        MT = sdr[lowest_indices, :]
        depth_GS = depth_GS[lowest_indices]
        M0 = M0_GS[lowest_indices]

        """ Calculate take-off angles"""
        takeoff_angles = ["P", "S", "pP"]
        angles = []
        for phase in takeoff_angles:
            angles.append(
                fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance, takeoffs=True)
            )
        """ Beachball plot """
        fig = Plot_GS_BB(
            MT[:, 0],
            MT[:, 1],
            MT[:, 2],
            azimuths=[event.az, event.az, event.az],
            inc_angles=angles,
            phase_names=takeoff_angles,
            color=color_plot,
        )
        plt.savefig(
            pjoin(
                h5_file_folder,
                f"GS_BBB__{event.name}_{depth}_{misfit_name}_{fwd.veloc_name}_Post.svg",
            ),
            dpi=300,
        )
        plt.close()

        """ Generate Green's functions per depth """
        syn_tts = []
        syn_GFs = []
        for i, phase in enumerate(phases):
            syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)

            syn_GF = fwd.get_greens_functions(
                comp=components[i],
                depth=depth,
                distance=event.distance,
                lat_src=event.latitude,
                lon_src=event.longitude,
                rec=rec,
                tstar=tstars[i],
                LQT=LQT_value,
                inc=inc,
                baz=baz,
                M0=1.0,
                filter=filter_par,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
            )
            syn_GFs.append(syn_GF)
            syn_tts.append(syn_tt)

        if plot_extra_phases is not None:
            extra_arrs = []
            for j, extraphase in enumerate(plot_extra_phases):
                arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
                extra_arrs.append(arr)
        else:
            extra_arrs = None

        fig, ax = waveform_plot(
            syn_GFs=syn_GFs,
            syn_tts=syn_tts,
            obs_tts=obs_tt,
            fwd=fwd,
            misfit_weight_len=misfit_weight_len,
            event=event,
            phases=phases,
            components=components,
            t_pre=t_pre,
            t_post=t_post,
            MTs=MT,
            M0s=M0,
            fmin=fmin,
            fmax=fmax,
            zerophase=zerophase,
            plot_extra_phases=plot_extra_phases,
            extra_arrs=extra_arrs,
            color_plot=color_plot,
            Ylims=Ylims,
        )

        """ DIRECT: """
        h5_file_path = pjoin(
            h5_file_folder,
            f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
        )

        (
            depth_Direct,
            MT,
            DC_MT,
            CLVD_MT,
            misfit_L2_Direct,
            Epsilon,
            M0,
            M0_DC,
            M0_CLVD,
            angles,
            cond_nr,
        ) = _ReadH5.Read_Direct_Inversion(h5_file_path, amount_of_phases=5)
        Total_L2_Direct = np.sum(misfit_L2_Direct)

        color_plot = "r"

        MT_FULL_ = MT  # np.array([MT[0], MT[2], MT[1], MT[4], MT[3], MT[5],])
        DC_MT_ = DC_MT  # np.array([DC_MT[0], DC_MT[2], DC_MT[1], DC_MT[4], DC_MT[3], DC_MT[5],])
        CLVD_MT_ = CLVD_MT  # np.array(
        #     [CLVD_MT[0], CLVD_MT[2], CLVD_MT[1], CLVD_MT[4], CLVD_MT[3], CLVD_MT[5],]
        # )

        fig = Plot_Direct_BB(
            MT_Full=MT_FULL_ / M0,
            Eps=Epsilon,
            MT_DC=DC_MT_ / M0_DC,
            M0_DC=M0_DC,
            MT_CLVD=CLVD_MT_ / M0_CLVD,
            M0_CLVD=M0_CLVD,
            azimuths=[event.az, event.az, event.az],
            inc_angles=angles,
            phase_names=["P", "S", "pP"],
            color=color_plot,
            height=19.0,
            horizontal=True,
        )

        plt.savefig(
            pjoin(
                h5_file_folder,
                f"Direct_BB_{event.name}_{depth}_{misfit_name}_{fwd.veloc_name}_Post.svg",
            ),
            dpi=300,
        )
        plt.close()

        MT = np.expand_dims(DC_MT, axis=0)
        M0 = np.expand_dims(M0_DC, axis=0)

        fig = waveform_plot(
            syn_GFs=syn_GFs,
            syn_tts=syn_tts,
            obs_tts=obs_tt,
            fwd=fwd,
            misfit_weight_len=misfit_weight_len,
            event=event,
            phases=phases,
            components=components,
            t_pre=t_pre,
            t_post=t_post,
            MTs=MT,
            M0s=M0,
            fmin=fmin,
            fmax=fmax,
            zerophase=zerophase,
            plot_extra_phases=plot_extra_phases,
            extra_arrs=extra_arrs,
            color_plot=color_plot,
            Ylims=Ylims,
            fig=fig,
            ax=ax,
        )

        if Return_Fig:
            return fig
        else:
            plt.savefig(
                pjoin(
                    h5_file_folder,
                    f"COMBINED_waveforms_{event.name}_{depth}_{misfit_name}_{fwd.veloc_name}_Post.svg",
                ),
                dpi=300,
            )
            plt.close()


def waveform_plot_copy(
    syn_GFs: obspy.Stream,
    syn_tts: [float],
    obs_tts: [float],
    fwd: _Forward._AbstractForward,
    misfit_weight_len: float,
    event: obspy.core.event.Event,
    phases: [str],
    components: [float],
    t_pre: [float],
    t_post: [float],
    MTs: _List[float],
    M0s: [float],
    depth: float,
    rec: instaseis.Receiver,
    tstar: _Union[float, str],
    phase_corrs: [float],
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = None,
    plot_extra_phases: [str] = None,
    extra_arrs: [float] = None,
    color_plot: str = None,
    Ylims: [float] = None,
    fig: [bool] = None,
    ax: [bool] = None,
):
    """ Waveform plot """

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    inv_phases = len(phases)

    """ Prepare all phases for the plot """
    all_phase_combs = set(["PZ", "PR", "SZ", "SR", "ST"])
    pc = set([p + c for p, c in zip(phases, components)])
    PorS = set([phase[0] for phase in pc])
    missing = all_phase_combs - pc
    phases_to_plot = []
    comps_to_plot = []
    corrs_missing = []
    tstar_missing = []
    t_pres_missing = []
    t_posts_missing = []
    Ylims_missing = []
    for missing_phase in missing:
        if missing_phase[0] in PorS:
            corrs_missing.append(phase_corrs[phases.index(missing_phase[0])])
            tstar_missing.append(tstar[phases.index(missing_phase[0])])
            t_pres_missing.append(t_pre[phases.index(missing_phase[0])])
            t_posts_missing.append(t_post[phases.index(missing_phase[0])])
            Ylims_missing.append(Ylims[phases.index(missing_phase[0])])
            phases_to_plot.append(missing_phase[0])
            comps_to_plot.append(missing_phase[1])
    """ Add missing observed arrival times"""
    for i, phase in enumerate(phases_to_plot):
        obs_tts.append(utct(event.picks[phase]) - event.origin_time + corrs_missing[i])
    """ Add missing synthetic GF """
    for i, phase in enumerate(phases_to_plot):
        syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)
        syn_GF = fwd.get_greens_functions(
            comp=comps_to_plot[i],
            depth=depth,
            distance=event.distance,
            lat_src=event.latitude,
            lon_src=event.longitude,
            rec=rec,
            tstar=tstar_missing[i],
            LQT=False,
            inc=None,
            baz=None,
            M0=1,
            filter=filter_par,
            fmin=fmin,
            fmax=fmax,
            zerophase=zerophase,
        )
        syn_GFs.append(syn_GF)
        syn_tts.append(syn_tt)

    phases_to_plot = phases + phases_to_plot
    comps_to_plot = components + comps_to_plot
    t_pres_missing = t_pre + t_pres_missing
    t_posts_missing = t_post + t_posts_missing
    Ylims_missing = Ylims + Ylims_missing

    if fig is None and ax is None:
        # fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(18, 16))
        fig, ax = plt.subplots(nrows=len(phases_to_plot), ncols=1, sharex="all", figsize=(18, 20))

    for i, phase in enumerate(phases_to_plot):
        for n in range(len(M0s)):
            tr_syn_full = fwd.generate_synthetic_data(
                st_GF=syn_GFs[i], focal_mech=MTs[n], M0=M0s[n], slice=False,
            )

            tr_slice = tr_syn_full.slice(
                starttime=fwd.or_time + syn_tts[i] - t_pres_missing[i],
                endtime=fwd.or_time + syn_tts[i] + t_posts_missing[i],
            )
            if i < inv_phases:
                if n == 0:
                    if color_plot == "blue":
                        ax[i].plot(
                            tr_slice.times() - t_pres_missing[i],
                            tr_slice.data,
                            lw=2,
                            c=color_plot,
                            label="Synthetic (Grid-search)",
                        )
                    else:
                        ax[i].plot(
                            tr_slice.times() - t_pres_missing[i],
                            tr_slice.data,
                            lw=3,
                            c=color_plot,
                            label="Synthetic (Direct)",
                        )

                else:
                    ax[i].plot(
                        tr_slice.times() - t_pres_missing[i], tr_slice.data, lw=2, c=color_plot,
                    )
            ax[i].plot(
                tr_syn_full.times() - (syn_tts[i] - fwd.start_cut),
                tr_syn_full.data,
                lw=1,
                c=color_plot,
            )
            # ax[i].legend()

            st = obspy.Stream()
            st += tr_slice

            if n == len(M0s) - 1:
                st_obs_full, sigmasPLOT = _PreProcess.prepare_event_data(
                    event=event,
                    phases=phases_to_plot,
                    components=comps_to_plot,
                    slice=False,
                    filter=filter_par,
                    fmin=fmin,
                    fmax=fmax,
                    zerophase=zerophase,
                    noise_level=False,
                )

                st_obs = st_obs_full.slice(
                    starttime=fwd.or_time + obs_tts[i] - t_pres_missing[i],
                    endtime=fwd.or_time + obs_tts[i] + t_posts_missing[i],
                )
                ax[i].plot(
                    st_obs_full[i].times() - obs_tts[i], st_obs_full[i].data, lw=1, c="k",
                )
                if i < inv_phases:
                    ax[i].plot(
                        st_obs[i].times() - t_pres_missing[i],
                        st_obs[i].data,
                        lw=3,
                        c="k",
                        label="Observed",
                    )

                st += st_obs[i]
                if Ylims_missing is None:
                    ax[i].set_ylim(global_min, global_max)
                else:
                    ax[i].set_ylim(-Ylims_missing[i], Ylims_missing[i])

                global_max = max([tr.data.max() for tr in st]) * 1.2
                global_min = min([tr.data.min() for tr in st]) * 1.2
                if i < inv_phases:
                    for axis in ["top", "bottom", "left", "right"]:
                        ax[i].spines[axis].set_linewidth(5)

                    ax[i].axvline(x=t_posts_missing[i], c="grey", ls="dashed")
                    ax[i].axvline(x=-t_pres_missing[i], c="grey", ls="dashed")
                    ax[i].axvspan(
                        -t_pres_missing[i], misfit_weight_len, facecolor="grey", alpha=0.3
                    )
                    ax[i].axvspan(
                        misfit_weight_len, t_posts_missing[i], facecolor="grey", alpha=0.1
                    )
                ymin = ax[i].get_ylim()[0]
                ymax = ax[i].get_ylim()[1]
                if event.name == "S0173a" and phase == "P":
                    ax[i].axvspan(
                        t_posts_missing[i], t_posts_missing[i] + 23, facecolor="green", alpha=0.1
                    )
                    ax[i].text(
                        29,
                        ymin * 0.8,
                        "Glitch",
                        verticalalignment="center",
                        color="green",
                        fontsize=35,
                    )

                ax[i].axvline(x=0.0, c="dimgrey", lw=2)
                ax[i].text(
                    0 + 0.1,
                    ymax * 0.8,
                    phase,
                    verticalalignment="center",
                    color="dimgray",
                    fontsize=30,
                )
                ax[i].text(
                    s="%s%s" % (phases_to_plot[i], comps_to_plot[i]),
                    x=0.99,
                    y=0.75,
                    ha="right",
                    transform=ax[i].transAxes,
                    color=color_plot,
                    fontsize=40,
                )

                # Extra phase arrivals:
                if plot_extra_phases is not None:
                    for j, extraphase in enumerate(plot_extra_phases):
                        arr = extra_arrs[j]
                        if arr:
                            if arr - syn_tts[i] > 0 and arr - syn_tts[i] < 31:
                                ax[i].axvline(x=arr - syn_tts[i], c="dimgrey", lw=2)
                                ax[i].text(
                                    arr - syn_tts[i] + 0.1,
                                    ymax * 0.75,
                                    extraphase,
                                    verticalalignment="center",
                                    color="dimgrey",
                                    fontsize=30,
                                    rotation=90,
                                )
                ax[i].tick_params(axis="both", which="major", labelsize=35)
                ax[i].tick_params(axis="both", which="minor", labelsize=25)

                ax[i].get_yaxis().get_offset_text().set_visible(False)
                ax_max = max(ax[i].get_yticks())
                exponent_axis = np.floor(np.log10(ax_max)).astype(int)
                # ax[i].annotate(
                #     r"$\times$10$^{%i}$" % (exponent_axis),
                #     xy=(0.01, 0.75),
                #     xycoords="axes fraction",
                #     fontsize=32,
                # )

    fig.text(
        0.9,
        0.88,
        "M0: %.2e" % (M0s[0]),
        ha="right",
        va="bottom",
        size="medium",
        color="black",
        fontsize=40,
    )
    fig.text(0.01, 0.5, "Displacement (nm)", va="center", rotation="vertical", fontsize=45)
    fig.text(
        0.5,
        0.88,
        event.name,
        ha="center",
        va="bottom",
        size="x-large",
        color=color_plot,
        fontsize=45,
    )

    ax[0].legend(
        prop={"size": 35},
        loc="center left",
        bbox_to_anchor=(0.12, 0.93),
        bbox_transform=fig.transFigure,
    )

    ax[-1].set_xlim(-10.0, 32.0)
    ax[-1].set_xlabel("time after phase (s)", fontsize=45)
    return fig, ax


def waveform_plot(
    syn_GFs: obspy.Stream,
    syn_tts: [float],
    obs_tts: [float],
    fwd: _Forward._AbstractForward,
    misfit_weight_len: float,
    event: obspy.core.event.Event,
    phases: [str],
    components: [float],
    t_pre: [float],
    t_post: [float],
    MTs: _List[float],
    M0s: [float],
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = None,
    plot_extra_phases: [str] = None,
    extra_arrs: [float] = None,
    color_plot: str = None,
    Ylims: [float] = None,
    fig: [bool] = None,
    ax: [bool] = None,
):
    """ Waveform plot """

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    if fig is None and ax is None:
        # fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(18, 16))
        fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(18, 20))

    for i, phase in enumerate(phases):
        for n in range(len(M0s)):
            tr_syn_full = fwd.generate_synthetic_data(
                st_GF=syn_GFs[i], focal_mech=MTs[n], M0=M0s[n], slice=False,
            )

            tr_slice = tr_syn_full.slice(
                starttime=fwd.or_time + syn_tts[i] - t_pre[i],
                endtime=fwd.or_time + syn_tts[i] + t_post[i],
            )

            if n == 0:
                if color_plot == "b":
                    ax[i].plot(
                        tr_slice.times() - t_pre[i],
                        tr_slice.data,
                        lw=2,
                        c=color_plot,
                        label="Synthetic (Grid-search)",
                    )
                else:
                    ax[i].plot(
                        tr_slice.times() - t_pre[i],
                        tr_slice.data,
                        lw=3,
                        c=color_plot,
                        label="Synthetic (Direct)",
                    )

            else:
                ax[i].plot(
                    tr_slice.times() - t_pre[i], tr_slice.data, lw=2, c=color_plot,
                )
            ax[i].plot(
                tr_syn_full.times() - (syn_tts[i] - fwd.start_cut),
                tr_syn_full.data,
                lw=1,
                c=color_plot,
            )
            # ax[i].legend()

            st = obspy.Stream()
            st += tr_slice

            if n == len(M0s) - 1:
                st_obs_full, sigmasPLOT = _PreProcess.prepare_event_data(
                    event=event,
                    phases=phases,
                    components=components,
                    slice=False,
                    filter=filter_par,
                    fmin=fmin,
                    fmax=fmax,
                    zerophase=zerophase,
                    noise_level=False,
                )

                st_obs = st_obs_full.slice(
                    starttime=fwd.or_time + obs_tts[i] - t_pre[i],
                    endtime=fwd.or_time + obs_tts[i] + t_post[i],
                )
                ax[i].plot(
                    st_obs_full[i].times() - obs_tts[i], st_obs_full[i].data, lw=1, c="k",
                )
                ax[i].plot(
                    st_obs[i].times() - t_pre[i], st_obs[i].data, lw=3, c="k", label="Observed",
                )

                st += st_obs[i]
                if Ylims is None:
                    ax[i].set_ylim(global_min, global_max)
                else:
                    ax[i].set_ylim(-Ylims[i], Ylims[i])

                global_max = max([tr.data.max() for tr in st]) * 1.2
                global_min = min([tr.data.min() for tr in st]) * 1.2
                ax[i].axvline(x=t_post[i], c="grey", ls="dashed")
                ax[i].axvline(x=-t_pre[i], c="grey", ls="dashed")
                ax[i].axvspan(-t_pre[i], misfit_weight_len, facecolor="grey", alpha=0.2)
                ymin = ax[i].get_ylim()[0]
                ymax = ax[i].get_ylim()[1]
                if event.name == "S0173a" and phase == "P":
                    ax[i].axvspan(t_post[i], t_post[i] + 23, facecolor="green", alpha=0.1)
                    ax[i].text(
                        29,
                        ymin * 0.8,
                        "Glitch",
                        verticalalignment="center",
                        color="green",
                        fontsize=35,
                    )

                ax[i].axvline(x=0.0, c="dimgrey", lw=2)
                ax[i].text(
                    0 + 0.1,
                    ymax * 0.8,
                    phase,
                    verticalalignment="center",
                    color="dimgray",
                    fontsize=30,
                )
                ax[i].text(
                    s="%s%s" % (phases[i], components[i]),
                    x=0.99,
                    y=0.75,
                    ha="right",
                    transform=ax[i].transAxes,
                    color=color_plot,
                    fontsize=40,
                )

                # Extra phase arrivals:
                if plot_extra_phases is not None:
                    for j, extraphase in enumerate(plot_extra_phases):
                        arr = extra_arrs[j]
                        if arr:
                            if arr - syn_tts[i] > 0 and arr - syn_tts[i] < 31:
                                ax[i].axvline(x=arr - syn_tts[i], c="dimgrey", lw=2)
                                ax[i].text(
                                    arr - syn_tts[i] + 0.1,
                                    ymax * 0.75,
                                    extraphase,
                                    verticalalignment="center",
                                    color="dimgrey",
                                    fontsize=30,
                                    rotation=90,
                                )
                ax[i].tick_params(axis="both", which="major", labelsize=35)
                ax[i].tick_params(axis="both", which="minor", labelsize=25)

                ax[i].get_yaxis().get_offset_text().set_visible(False)
                ax_max = max(ax[i].get_yticks())
                exponent_axis = np.floor(np.log10(ax_max)).astype(int)
                # ax[i].annotate(
                #     r"$\times$10$^{%i}$" % (exponent_axis),
                #     xy=(0.01, 0.75),
                #     xycoords="axes fraction",
                #     fontsize=32,
                # )

    fig.text(
        0.9,
        0.88,
        "M0: %.2e" % (M0s[0]),
        ha="right",
        va="bottom",
        size="medium",
        color="black",
        fontsize=40,
    )
    fig.text(0.01, 0.5, "Displacement (nm)", va="center", rotation="vertical", fontsize=45)
    fig.text(
        0.5,
        0.88,
        event.name,
        ha="center",
        va="bottom",
        size="x-large",
        color=color_plot,
        fontsize=45,
    )

    ax[0].legend(
        prop={"size": 35},
        loc="center left",
        bbox_to_anchor=(0.12, 0.93),
        bbox_transform=fig.transFigure,
    )

    ax[-1].set_xlim(-10.0, 32.0)
    ax[-1].set_xlabel("time after phase (s)", fontsize=45)
    return fig, ax


def Source_Uncertainty_OLD(
    h5_file_folder: str,
    event_name: str,
    method: str,
    misfit_name: str,
    fwd: _Forward._AbstractForward,
    phases: [str],
    components: [str],
    depths: [float],
    DOF: float,
    fmin: float = None,
    fmax: float = None,
):

    for idepth, depth in enumerate(depths):
        print(depth)
        # if method == "GS":
        color_plot = "b"
        h5_file_path = pjoin(
            h5_file_folder,
            f"GS_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
        )
        depth_GS, sdr, M0_GS, misfit_L2_GS = _ReadH5.Read_GS_h5(
            Filename=h5_file_path, amount_of_phases=len(phases)
        )
        Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
        n_lowest = 50
        # n_lowest = int(len(Total_L2_GS) * 0.05)
        # lowest_indices = Total_L2_GS.argsort()[0:n_lowest:50]
        lowest_indices = Total_L2_GS.argsort()[0:n_lowest]

        GOF_GS = (Total_L2_GS / DOF)[lowest_indices]
        M0 = M0_GS[lowest_indices]

        sdrs = sdr[lowest_indices, :]
        MT_Full = np.zeros((sdrs.shape[0], 6))
        for i in range(MT_Full.shape[0]):
            MT_Full[i, :] = _GreensFunctions.convert_SDR(sdrs[i, 0], sdrs[i, 1], sdrs[i, 2], M0[i])
            MT_Full[i, 3] = -MT_Full[i, 4]
            MT_Full[i, 5] = -MT_Full[i, 5]

        if idepth == 0:
            M0_plot_GS = M0
            MT_GS = MT_Full
            MT_sdrs = sdrs
            weights_GS = np.exp(-GOF_GS)
        else:
            M0_plot_GS = np.hstack((M0_plot_GS, M0))
            MT_GS = np.vstack((MT_GS, MT_Full))
            MT_sdrs = np.vstack((MT_sdrs, sdrs))
            weights_GS = np.hstack((weights_GS, np.exp(-GOF_GS)))
        (values, counts) = np.unique(sdrs[:, 0], return_counts=True)
        ind = np.argmax(counts)
        print("Strike:", values[ind])
        (values, counts) = np.unique(sdrs[:, 1], return_counts=True)
        ind = np.argmax(counts)
        print("Dip:", values[ind])
        (values, counts) = np.unique(sdrs[:, 2], return_counts=True)
        ind = np.argmax(counts)
        print("Rake:", values[ind])

        # else:
        if event_name == "S0183a":
            pass
        else:
            h5_file_path = pjoin(
                h5_file_folder,
                f"Direct_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )

            (
                depth_Direct,
                MT_Full,
                DC_MT,
                CLVD_MT,
                misfit_L2_Direct,
                Epsilon,
                M0_Direct,
                M0_DC,
                M0_CLVD,
                angles,
                cond_nr,
            ) = _ReadH5.Read_Direct_Inversion(h5_file_path, amount_of_phases=len(phases))
            Total_L2_Direct = np.sum(misfit_L2_Direct)
            GOF_Direct = Total_L2_Direct / DOF
            DC_MT = np.expand_dims(DC_MT, axis=0)
            # DC_MT = DC_MT / M0_DC
            DC_MT[0, 3] = -DC_MT[0, 4]
            DC_MT[0, 5] = -DC_MT[0, 5]
            if idepth == 0:
                M0_plot_Direct = M0_DC
                MT_Direct = DC_MT
                weights_Direct = np.exp(-GOF_Direct)
            else:
                M0_plot_Direct = np.hstack((M0_plot_Direct, M0_DC))
                MT_Direct = np.vstack((MT_Direct, DC_MT))
                weights_Direct = np.hstack((weights_Direct, np.exp(-GOF_Direct)))

            color_plot = "r"

    # MT_names = ["mrr", "mtt", "mpp", "mrt", "mrp", "mtp"]
    MT_names = ["mzz", "mxx", "myy", "mxz", "myz", "mxy"]

    ## =====================BIG PLOT ===================
    fig, ax = plt.subplots(
        nrows=2, ncols=6, figsize=(12, 6), sharey="row", gridspec_kw={"height_ratios": [3, 1]}
    )
    # for i in range(6):
    #     hist_1 = MT_GS[:, i] / np.max(MT_GS[:, i] )
    #     ax[0, i].hist(
    #         hist_1,
    #         bins=18,
    #         alpha=0.4,
    #         color="steelblue",
    #         edgecolor="none",
    #         label="GS",
    #         weights=weights_GS,
    #         density=True,
    #     )
    #     hist_2 = MT_Direct[:, i] / np.max(MT_Direct[:, i] )
    #     ax[0, i].hist(
    #         hist_2,
    #         bins=18,
    #         alpha=0.4,
    #         color="red",
    #         edgecolor="none",
    #         label="Direct",
    #         weights=weights_Direct,
    #         density=True,
    #     )
    #     ax[0, i].set_xlabel(MT_names[i], fontsize=18)
    #     ax[0, i].tick_params(axis="x", labelsize=15)
    #     ax[0, i].tick_params(axis="y", labelsize=15)
    #     ax[0, i].ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
    #     ax[0, i].set_xlim(-1, 1)

    #     h = np.histogram(MT_GS[:, i],bins=18)
    #     h = np.vstack((0.5*(h[1][:-1]+h[1][1:]),h[0])).T
    #     array_convt = h
    #     km = KMeans(n_clusters=1).fit(array_convt)
    #     clusters = km.predict(array_convt)
    #     kmeans = km.cluster_centers_
    #     c_mean_distances = []
    #     for i, (cx,cy) in enumerate(kmeans):
    #             mean_distance = k_mean_distance_2d(array_convt, cx,cy,i, clusters)
    #             c_mean_distances.append(mean_distance)

    #     ax[0, i].axvline(x=kmeans[0][0] + c_mean_distances[0], c="blue", lw=1, ls="--")
    #     ax[0, i].axvline(x=kmeans[0][0] - c_mean_distances[0], c="blue", lw=1, ls="--")
    #     ax[0, i].axvline(x=kmeans[0][0] - c_mean_distances[0], c="blue", lw=1)

    #     n_GS, bins_GS = np.histogram(MT_GS[:, i], bins=18, weights=weights_GS, density=True)
    #     mids = 0.5 * (bins_GS[1:] + bins_GS[:-1])
    #     mean_GS = np.average(mids, weights=n_GS)
    #     # (values, counts) = np.unique(MT_GS[:, i], return_counts=True)
    #     # ind = np.argmax(counts)
    #     # mean_GS = values[ind]
    #     var_GS = np.sqrt(np.average((mids - mean_GS) ** 2, weights=n_GS))
    #     # ax[0, i].axvline(x=mean_GS, c="steelblue", lw=1)
    #     # ax[0, i].axvline(x=mean_GS + var_GS, c="steelblue", lw=1, ls="--")
    #     # ax[0, i].axvline(x=mean_GS - var_GS, c="steelblue", lw=1, ls="--")

    #     n_Direct, bins_Direct = np.histogram(
    #         MT_Direct[:, i], bins=18, weights=weights_Direct, density=True
    #     )
    #     mids = 0.5 * (bins_Direct[1:] + bins_Direct[:-1])
    #     mean_Direct = np.average(mids, weights=n_Direct)
    #     var_Direct = np.sqrt(np.average((mids - mean_Direct) ** 2, weights=n_Direct))
    #     # ax[0, i].axvline(x=mean_Direct, c="red", lw=1)
    #     # ax[0, i].axvline(x=mean_Direct + var_Direct, c="red", lw=1, ls="--")
    #     # ax[0, i].axvline(x=mean_Direct - var_Direct, c="red", lw=1, ls="--")

    # ax[0, 0].set_ylabel("Frequency", fontsize=18)
    # # ax[0,0].set_ylim(0, 10)
    # ax[0, 0].legend()

    sdr_names = ["strike", "dip", "rake"]
    sdr_mins = [0, 0, -180]
    sdr_maxs = [360, 90, 180]
    ## ========================================

    # strike_rad = MT_sdrs[:, 0]
    # dip_rad = MT_sdrs[:, 1]
    # n = np.array(
    #     [
    #         -np.sin(dip_rad) * np.sin(strike_rad),
    #         -np.sin(dip_rad) * np.cos(strike_rad),
    #         np.cos(dip_rad),
    #     ]
    # )
    # x = n[0, :]
    # y = n[1, :]
    # z = n[2, :]

    # km_xyz = KMeans(n_clusters=2).fit(n.T)
    # clusters_xyz = km_xyz.predict(n.T)
    # kmeans_xyz = km_xyz.cluster_centers_
    # c_mean_distances_xyz = []
    # for i, (cx, cy, cz) in enumerate(kmeans_xyz):
    #     mean_distance = k_mean_distance_3d(n.T, cx, cy, cz, i, clusters_xyz)
    #     c_mean_distances_xyz.append(mean_distance)

    # dips = []
    # strikes = []
    # for cluster in kmeans_xyz:
    #     s = np.rad2deg(np.arctan2(cluster[1], cluster[0]))
    #     wrap = s % 360
    #     st = 360.0 - 90.0 - wrap
    #     s_new = st % 360.0

    #     dip = np.rad2deg(np.arctan2(np.sqrt(cluster[0] ** 2 + cluster[1] ** 2), cluster[2]))
    #     if dip >= 90 or dip < 0:
    #         if dip == 90:
    #             dip = 89.0
    #             strike = (s_new + 180.0) % 360.0
    #             # if rake == -180:
    #             #     pass
    #             # else:
    #             #     rake = -rake
    #         else:
    #             d_new = 90.0 - dip
    #             dip = d_new % 90.0
    #             strike = (s_new + 180.0) % 360.0
    #             # if rake == -180:
    #             #     pass
    #             # else:
    #             #     rake = -rake

    #     dips.append(dip)
    #     strikes.append(strike)

    # h = np.histogram(MT_sdrs[:, 0],bins=18)
    # h = np.vstack((0.5*(h[1][:-1]+h[1][1:]),h[0])).T
    # array_convt = h
    # strike_rad = np.deg2rad(MT_sdrs[:, 0])
    # strike_unwrap = np.unwrap(strike_rad)

    # array_convt=strike_rad.reshape(len(strike_rad),1)
    array_convt = MT_sdrs[:, 0].reshape(len(MT_sdrs[:, 0]), 1)
    km_s = KMeans(n_clusters=2).fit(array_convt)
    clusters_s = km_s.predict(array_convt)
    kmeans_s = km_s.cluster_centers_
    c_mean_distances_s = []
    # for i, (cx,cy) in enumerate(kmeans_s):
    #         mean_distance = k_mean_distance_2d(h, cx,cy,i, clusters_s)
    #         c_mean_distances_s.append(mean_distance)
    for i, cx in enumerate(kmeans_s):
        mean_distance = k_mean_distance(array_convt, cx, i, clusters_s)
        c_mean_distances_s.append(mean_distance)

    array_convt = MT_sdrs[:, 1].reshape(len(MT_sdrs[:, 1]), 1)
    km_d = KMeans(n_clusters=2).fit(array_convt)
    clusters_d = km_d.predict(array_convt)
    kmeans_d = km_d.cluster_centers_
    c_mean_distances_d = []
    for i, cx in enumerate(kmeans_d):
        mean_distance = k_mean_distance(array_convt, cx, i, clusters_d)
        c_mean_distances_d.append(mean_distance)

    array_convt = MT_sdrs[:, 2]
    array_convt = array_convt[array_convt < 0]
    array_convt = array_convt.reshape(len(array_convt), 1)
    km_r = KMeans(n_clusters=2).fit(array_convt)
    clusters_r = km_r.predict(array_convt)
    kmeans_r = km_r.cluster_centers_
    c_mean_distances_r = []
    for i, cx in enumerate(kmeans_r):
        mean_distance = k_mean_distance(array_convt, cx, i, clusters_r)
        c_mean_distances_r.append(mean_distance)

    c_mean_distances = [c_mean_distances_s, c_mean_distances_d, c_mean_distances_r]
    print(c_mean_distances)

    strike_aux, dip_aux, rake_aux = aux_plane(kmeans_s[0][0], kmeans_d[0][0], kmeans_r[0][0])
    # strike_aux, dip_aux, rake_aux = aux_plane(strikes[0], dips[0], kmeans_r[0][0])
    f_planes = [kmeans_s[0][0], kmeans_d[0][0], kmeans_r[0][0]]
    # f_planes = [strikes[0], dips[0], kmeans_r[0][0]]
    kmeans_tot = [kmeans_s[1][0], kmeans_d[1][0], kmeans_r[1][0]]
    # kmeans_tot = [strikes[1],dips[1], kmeans_r[1][0]]
    print(f"plane 1: {f_planes}")
    aux_planes = [strike_aux, dip_aux, rake_aux]
    print(f"plane 2: {aux_planes}")

    # MT_planes = (
    #     _GreensFunctions.convert_SDR(kmeans_s[0][0], kmeans_d[0][0], kmeans_r[0][0], M0_DC) / M0_DC
    # )
    # MT_planes[3] = -MT_planes[3]
    # MT_planes[5] = -MT_planes[5]
    ## =====================BIG PLOT ===================
    # for i in range(6):
    #     ax[0, i].axvline(x=MT_planes[i], c="red", lw=1, label="Fault planes")

    # for axs in ax[1, :]:
    #     axs.remove()
    # for i in range(3):
    #     plot_nr = i * 2

    #     gs = ax[1, plot_nr].get_gridspec()
    #     # remove the underlying axes

    #     axbig = fig.add_subplot(gs[1, plot_nr : plot_nr + 2])
    #     # axbig.annotate('Big Axes \nGridSpec[1:, -1]', (0.1, 0.5)
    #     axbig.hist(
    #         MT_sdrs[:, i],
    #         bins=18,
    #         alpha=0.4,
    #         color="steelblue",
    #         edgecolor="none",
    #         weights=weights_GS,
    #         density=True,
    #     )
    #     axbig.axvline(x=f_planes[i], c="blue", lw=1, label="K-means + fault plane")
    #     axbig.axvline(x=kmeans_tot[i], c="green", lw=1, label="K-means")
    #     axbig.axvline(x=aux_planes[i], c="red", lw=1, label="Aux plane")

    #     axbig.axvline(x=f_planes[i] + c_mean_distances[i][0], c="blue", lw=1, ls="--")
    #     axbig.axvline(x=f_planes[i] - c_mean_distances[i][0], c="blue", lw=1, ls="--")
    #     axbig.axvline(x=kmeans_tot[i] + c_mean_distances[i][1], c="green", lw=1, ls="--")
    #     axbig.axvline(x=kmeans_tot[i] - c_mean_distances[i][1], c="green", lw=1, ls="--")

    #     axbig.set_xlabel(sdr_names[i], fontsize=18)
    #     axbig.tick_params(axis="x", labelsize=15)
    #     axbig.tick_params(axis="y", labelsize=15)
    #     axbig.ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
    #     axbig.set_xlim(sdr_mins[i], sdr_maxs[i])
    #     if i == 2:
    #         axbig.legend()
    fig1, ax1 = plt.subplots(nrows=1, ncols=3, figsize=(12, 3), sharey="row")
    for i in range(3):
        ax1[i].hist(
            MT_sdrs[:, i],
            bins=18,
            alpha=0.4,
            color="steelblue",
            edgecolor="none",
            weights=weights_GS,
            density=True,
        )
        # ax1[i].axvline(x=f_planes[i], c="blue", lw=1, label="fault plane 1")
        # ax1[i].axvline(x=kmeans_tot[i], c="green", lw=1, label="fault plane 2")
        # # ax1[i].axvline(x=aux_planes[i], c="red", lw=1, label = "Aux plane")

        # ax1[i].axvline(x=f_planes[i] + c_mean_distances[i][0], c="blue", lw=1, ls="--")
        # ax1[i].axvline(x=f_planes[i] - c_mean_distances[i][0], c="blue", lw=1, ls="--")
        # ax1[i].axvline(x=kmeans_tot[i] + c_mean_distances[i][1], c="green", lw=1, ls="--")
        # ax1[i].axvline(x=kmeans_tot[i] - c_mean_distances[i][1], c="green", lw=1, ls="--")

        ax1[i].set_xlabel(sdr_names[i], fontsize=18)
        ax1[i].tick_params(axis="x", labelsize=15)
        ax1[i].tick_params(axis="y", labelsize=15)
        ax1[i].ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
        ax1[i].set_xlim(sdr_mins[i], sdr_maxs[i])
        # if i == 2:
        #     ax1[i].legend()
    fig1.suptitle(event_name, fontsize=20)

    fig2, ax2 = plt.subplots(nrows=1, ncols=6, figsize=(12, 3), sharey="row")
    # MT_names = ["mrr", "mpp", "mtt", "mrp", "mrt", "mtp"]
    MT_names = ["mzz", "myy", "mxx", "myz", "mxz", "mxy"]

    # def integrator(f, data, freq):
    #     diffs = np.roll(data, -1) - data
    #     return (f(data[:-1]) * freq[:-1] * diffs[:-1]).sum()

    n, bins = np.histogram(M0_plot_GS, bins=18)  # , weights=weights_GS)
    mids = 0.5 * (bins[1:] + bins[:-1])

    mean = np.average(mids, weights=n)
    std = np.sqrt(np.average((mids - mean) ** 2, weights=n))
    print("M0", mean, std)
    print("MW", 2.0 / 3.0 * (np.log10(mean) - 9.1), 2.0 / 3.0 * (np.log10(std) - 9.1))

    for i in range(6):
        hist_1 = MT_GS[:, i]  # / np.max(MT_GS[:, i] )

        # freq_norm = weights_GS / integrator(lambda x: 1, hist_1, weights_GS)

        # mean = integrator(lambda x: x, hist_1, freq_norm)
        # std = integrator(lambda x: x ** 2, hist_1, freq_norm)

        # ax2[i].axvline(x=mean, c="steelblue", lw=1)
        # ax2[i].axvline(x=mean + std, c="steelblue", lw=1, ls="--")
        # ax2[i].axvline(x=mean - std, c="steelblue", lw=1, ls="--")

        n, bins = np.histogram(hist_1, bins=18, weights=weights_GS)
        mids = 0.5 * (bins[1:] + bins[:-1])

        mean = np.average(mids, weights=n)
        std = np.sqrt(np.average((mids - mean) ** 2, weights=n))
        ax2[i].axvline(x=mean, c="steelblue", lw=1, label="mean", alpha=0.5)
        ax2[i].axvline(x=mean + std, c="steelblue", lw=1, ls="--", label="std", alpha=0.5)
        ax2[i].axvline(x=mean - std, c="steelblue", lw=1, ls="--", alpha=0.5)

        print(MT_names[i], "%.2e" % mean, "%.2e" % std)
        # print(MT_names[i], hist_1.mean(), hist_1.std())
        ax2[i].hist(
            hist_1,
            bins=18,
            alpha=0.6,
            color="steelblue",
            edgecolor="none",
            label="GS",
            weights=weights_GS,
            density=True,
        )

        if event_name == "S0183a":
            pass
        else:
            hist_2 = MT_Direct[:, i]  # / np.max(MT_Direct[:, i] )
            # ax3 = ax2[i].twinx()
            ax2[i].hist(
                hist_2,
                bins=18,
                alpha=0.6,
                color="red",
                edgecolor="none",
                label="Direct",
                weights=weights_Direct,
                density=True,
            )
        # ax3.tick_params(axis="y", labelsize=15)
        # ax3.ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
        # ax3.tick_params(axis="y", labelcolor="red")
        # if i < 5:
        #     ax3.set_yticklabels([])
        # else:
        #     ax3.legend()

        ax2[i].set_xlabel(MT_names[i], fontsize=18)
        ax2[i].tick_params(axis="x", labelsize=15)
        ax2[i].tick_params(axis="y", labelsize=15)
        ax2[i].ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
        # ax2[i].tick_params(axis="y", labelcolor="steelblue")
        # ax2[i].set_xlim(-1, 1)
        # ax2[i].set_ylim(0, 0.05e-12)
        if event_name == "S0235b":
            ax2[i].set_ylim(0, 0.6e-12)
        elif event_name == "S0173a":
            ax2[i].set_ylim(0, 0.4e-12)
        elif event_name == "S0183a":
            ax2[i].set_ylim(0, 0.05e-12)

        ax2[i].get_xaxis().get_offset_text().set_visible(False)
        ax_max = max(ax2[i].get_xticks())
        if ax_max < 1.0:
            print("hello")
            ax_max = np.abs(ax2[i].get_xticks()[-2])
            print(ax_max)
        exponent_axis = np.floor(np.log10(ax_max)).astype(int)
        ax2[i].set_xlabel(
            f"{MT_names[i]}\n" + r"$\times$10$^{%i}$ (Nm)" % (exponent_axis), fontsize=14
        )
        # ax2[i].annotate(
        #     r"$\times$10$^{%i}$" % (exponent_axis), xy=(0.01, 0.82), xycoords="axes fraction",
        # )

        # h = np.histogram(MT_GS[:, i], bins=18)
        # h = np.vstack((0.5 * (h[1][:-1] + h[1][1:]), h[0])).T
        # array_convt = h
        # km = KMeans(n_clusters=1).fit(array_convt)
        # clusters = km.predict(array_convt)
        # kmeans = km.cluster_centers_
        # c_mean_distances = []
        # for j, (cx, cy) in enumerate(kmeans):
        #     mean_distance = k_mean_distance_2d(array_convt, cx, cy, i, clusters)
        #     c_mean_distances.append(mean_distance)

        # ax2[i].axvline(x=kmeans[0][0] + c_mean_distances[0], c="blue", lw=1, ls="--")
        # ax2[i].axvline(x=kmeans[0][0] - c_mean_distances[0], c="blue", lw=1, ls="--")
        # ax2[i].axvline(x=kmeans[0][0] - c_mean_distances[0], c="blue", lw=1)
    ax2[0].set_ylabel("Frequency", fontsize=18)
    # ax[0].set_ylim(0, 10)
    ax2[-1].legend()
    fig2.suptitle(event_name, fontsize=20)
    return fig2, fig1


def Source_Uncertainty(
    h5_file_folder: str,
    event_name: str,
    method: str,
    misfit_name: str,
    fwd: _Forward._AbstractForward,
    phases: [str],
    components: [str],
    depths: [float],
    DOF: float,
    fmin: float = None,
    fmax: float = None,
):
    for idepth, depth in enumerate(depths):
        print(depth)
        # if method == "GS":
        # else:
        if event_name == "S0183a":
            pass
        else:
            h5_file_path = glob.glob(
                pjoin(
                    h5_file_folder,
                    f"Direct_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}_*.hdf5",
                )
            )[0]

            (
                depth_Direct,
                MT_Full,
                DC_MT,
                CLVD_MT,
                misfit_L2_Direct,
                Epsilon,
                M0_Direct,
                M0_DC,
                M0_CLVD,
                angles,
                cond_nr,
            ) = _ReadH5.Read_Direct_Inversion(h5_file_path, amount_of_phases=5)
            if Epsilon > 0.2:
                continue
            if event_name == "S0173a" and depth > 60.0:
                continue
            print(f"chosen depth: {depth}")
            Total_L2_Direct = np.sum(misfit_L2_Direct)
            GOF_Direct = Total_L2_Direct / DOF
            DC_MT = np.expand_dims(DC_MT, axis=0)

            m = mtm.MomentTensor(
                mnn=DC_MT[0, 1],
                mee=DC_MT[0, 2],
                mdd=DC_MT[0, 0],
                mne=-DC_MT[0, 5],
                mnd=DC_MT[0, 3],
                med=-DC_MT[0, 4],
            )

            (s1, d1, r1), (s2, d2, r2) = m.both_strike_dip_rake()
            sdr_1 = np.expand_dims([s1, d1, r1], axis=0)
            sdr_2 = np.expand_dims([s2, d2, r2], axis=0)
            if "M0_plot_Direct" not in locals():
                M0_plot_Direct = M0_DC
                MT_Direct = DC_MT
                sdr_Direct1 = sdr_1
                sdr_Direct1 = np.vstack((sdr_Direct1, sdr_2))
                # sdr_Direct2 = sdr_2
                weights_Direct = np.exp(-GOF_Direct)
                weights_Direct = np.hstack((weights_Direct, np.exp(-GOF_Direct)))
                weights_Direct_MT = np.exp(-GOF_Direct)
            else:
                M0_plot_Direct = np.hstack((M0_plot_Direct, M0_DC))
                MT_Direct = np.vstack((MT_Direct, DC_MT))
                sdr_Direct1 = np.vstack((sdr_Direct1, sdr_1))
                sdr_Direct1 = np.vstack((sdr_Direct1, sdr_2))
                # sdr_Direct2 = np.vstack((sdr_Direct2, sdr_2))
                weights_Direct = np.hstack((weights_Direct, np.exp(-GOF_Direct)))
                weights_Direct = np.hstack((weights_Direct, np.exp(-GOF_Direct)))
                weights_Direct_MT = np.hstack((weights_Direct_MT, np.exp(-GOF_Direct)))

            color_plot = "r"
        color_plot = "b"
        h5_file_path = glob.glob(
            pjoin(
                h5_file_folder,
                f"GS_{event_name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}_*.hdf5",
            )
        )[0]
        depth_GS, sdr, M0_GS, misfit_L2_GS = _ReadH5.Read_GS_h5(
            Filename=h5_file_path, amount_of_phases=len(phases)
        )
        Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
        lowest_ind = Total_L2_GS.argsort()
        Total_L2_GS.sort()
        misfit_low = Total_L2_GS[:] - Total_L2_GS[0]
        uncert = 0.05 * Total_L2_GS[0]
        inds = np.where(misfit_low < uncert)
        lowest_indices = lowest_ind[inds][:10]
        GOF_GS = (np.sum(misfit_L2_GS, axis=1) / DOF)[lowest_indices]
        M0 = M0_GS[lowest_indices]

        sdrs = sdr[lowest_indices, :]
        MT_Full = np.zeros((sdrs.shape[0], 6))
        for i in range(MT_Full.shape[0]):
            MT_Full[i, :] = _GreensFunctions.convert_SDR(sdrs[i, 0], sdrs[i, 1], sdrs[i, 2], M0[i])

        if "M0_plot_GS" not in locals():
            M0_plot_GS = M0[:1]
            MT_GS = MT_Full
            MT_sdrs = sdrs
            weights_GS = np.exp(-GOF_GS)
        else:
            M0_plot_GS = np.hstack((M0_plot_GS, M0[:1]))
            MT_GS = np.vstack((MT_GS, MT_Full))
            MT_sdrs = np.vstack((MT_sdrs, sdrs))
            weights_GS = np.hstack((weights_GS, np.exp(-GOF_GS)))

    mean_Full = []
    std_Full = []
    n, bins = np.histogram(M0_plot_GS, bins=18)  # , weights=weights_GS)
    mids = 0.5 * (bins[1:] + bins[:-1])

    mean = np.average(mids, weights=n)
    std = np.sqrt(np.average((mids - mean) ** 2, weights=n))
    print("M0", mean, std)
    print("MW", 2.0 / 3.0 * (np.log10(mean) - 9.1), 2.0 / 3.0 * (np.log10(std) - 9.1))
    # from scipy import stats
    MT_names = ["mzz", "myy", "mxx", "mxz", "myz", "mxy"]
    for i in range(6):
        hist_1 = MT_GS[:, i]
        n, bins = np.histogram(hist_1, bins=18, weights=weights_GS)
        # hist_1 = MT_Direct[:, i]
        # n, bins = np.histogram(hist_1, bins=18, weights=weights_Direct_MT)

        mids = 0.5 * (bins[1:] + bins[:-1])

        mean = np.average(mids, weights=n)
        std = np.sqrt(np.average((mids - mean) ** 2, weights=n))
        print(MT_names[i], "%.2e" % mean, "%.2e" % std)
        # (values, counts) = np.unique(hist_1, return_counts=True)
        # ind = np.argmax(counts)
        # mean = values[ind]

        # n, bins = np.histogram(hist_1, bins=18, weights=weights_GS)
        # mids = 0.5 * (bins[1:] + bins[:-1])
        # mean = stats.mode(hist_1).mode[0]
        # if mean == -0.0:
        #     print("hallo")
        #     mean = 0.0
        # print(f"{mean}")
        # mean = np.average(hist_1, weights=weights_GS)
        # mean = np.mean(hist_1)
        mean_Full.append(mean)
        std_Full.append(np.sqrt(np.average((mids - mean) ** 2, weights=n)))

    m = mtm.MomentTensor(
        mnn=mean_Full[1],
        mee=mean_Full[2],
        mdd=mean_Full[0],
        mne=-mean_Full[5],
        mnd=mean_Full[3],
        med=-mean_Full[4],
    )

    (s1, d1, r1), (s2, d2, r2) = m.both_strike_dip_rake()
    mean_Full1 = [s1, d1, r1]
    print(mean_Full1)
    mean_Full2 = [s2, d2, r2]
    print(mean_Full2)

    m = mtm.MomentTensor(
        mnn=std_Full[1],
        mee=std_Full[2],
        mdd=std_Full[0],
        mne=-std_Full[5],
        mnd=std_Full[3],
        med=-std_Full[4],
    )

    (s1, d1, r1), (s2, d2, r2) = m.both_strike_dip_rake()
    std_Full1 = [s1, d1, r1]
    std_Full2 = [s2, d2, r2]

    """ K-means """
    array_convt = MT_sdrs[:, 0].reshape(len(MT_sdrs[:, 0]), 1)
    km_s = KMeans(n_clusters=2).fit(array_convt)
    clusters_s = km_s.predict(array_convt)
    kmeans_s = km_s.cluster_centers_
    c_mean_distances_s = []
    # for i, (cx,cy) in enumerate(kmeans_s):
    #         mean_distance = k_mean_distance_2d(h, cx,cy,i, clusters_s)
    #         c_mean_distances_s.append(mean_distance)
    for i, cx in enumerate(kmeans_s):
        mean_distance = k_mean_distance(array_convt, cx, i, clusters_s)
        c_mean_distances_s.append(mean_distance)

    array_convt = MT_sdrs[:, 1].reshape(len(MT_sdrs[:, 1]), 1)
    km_d = KMeans(n_clusters=2).fit(array_convt)
    clusters_d = km_d.predict(array_convt)
    kmeans_d = km_d.cluster_centers_
    c_mean_distances_d = []
    for i, cx in enumerate(kmeans_d):
        mean_distance = k_mean_distance(array_convt, cx, i, clusters_d)
        c_mean_distances_d.append(mean_distance)

    array_convt = MT_sdrs[:, 2]
    array_convt = array_convt[array_convt < 0]
    array_convt = array_convt.reshape(len(array_convt), 1)
    km_r = KMeans(n_clusters=2).fit(array_convt)
    clusters_r = km_r.predict(array_convt)
    kmeans_r = km_r.cluster_centers_
    c_mean_distances_r = []
    for i, cx in enumerate(kmeans_r):
        mean_distance = k_mean_distance(array_convt, cx, i, clusters_r)
        c_mean_distances_r.append(mean_distance)

    c_mean_distances = [c_mean_distances_s, c_mean_distances_d, c_mean_distances_r]
    print(c_mean_distances)

    strike_aux, dip_aux, rake_aux = aux_plane(kmeans_s[0][0], kmeans_d[0][0], kmeans_r[0][0])
    # strike_aux, dip_aux, rake_aux = aux_plane(strikes[0], dips[0], kmeans_r[0][0])
    f_planes = [kmeans_s[0][0], kmeans_d[0][0], kmeans_r[0][0]]
    # f_planes = [strikes[0], dips[0], kmeans_r[0][0]]
    kmeans_tot = [kmeans_s[1][0], kmeans_d[1][0], kmeans_r[1][0]]
    # kmeans_tot = [strikes[1],dips[1], kmeans_r[1][0]]
    print(f"plane 1: {f_planes}")
    aux_planes = [strike_aux, dip_aux, rake_aux]
    print(f"plane 2: {aux_planes}")

    """ """

    sdr_names = ["strike", "dip", "rake"]
    sdr_mins = [0, 0, -180]
    sdr_maxs = [360, 90, 180]
    ticklabels = [[0, 90, 180, 270, 360], [0, 15, 30, 45, 60, 75, 90], [-180, -90, 0, 90, 180]]

    # fig1, ax1 = plt.subplots(nrows=1, ncols=3, figsize=(12, 3), sharey="row")
    fig1, ax1 = plt.subplots(nrows=1, ncols=3, figsize=(28, 5))
    for i in range(3):
        ax1[i].hist(
            MT_sdrs[:, i],
            bins=18,
            alpha=0.8,
            range=(sdr_mins[i], sdr_maxs[i]),
            color="blue",
            edgecolor="none",
            weights=weights_GS,
            density=True,
            label="GS",
        )
        if event_name == "S0183a":
            pass
        else:
            ax2 = ax1[i].twinx()
            ax2.hist(
                sdr_Direct1[:, i],
                bins=18,
                range=(sdr_mins[i], sdr_maxs[i]),
                alpha=0.8,
                color="red",
                edgecolor="none",
                weights=weights_Direct,
                density=True,
                label="Direct",
            )
        # ax2.hist(
        #     sdr_Direct2[:, i],
        #     bins=18,
        #     alpha=0.4,
        #     color="darkred",
        #     edgecolor="none",
        #     weights=weights_Direct,
        #     density=True,
        # )
        ax1[i].axvline(x=mean_Full1[i], c="black", ls="dashed", lw=3, label="Fault plane 1")
        # ax1[i].axvline(
        #     x=mean_Full1[i] + std_Full1[i], c="red", lw=1, ls="--", label="std 1", alpha=0.5
        # )
        # ax1[i].axvline(x=mean_Full1[i] - std_Full1[i], c="red", lw=1, ls="--", alpha=0.5)

        ax1[i].axvline(x=mean_Full2[i], c="black", ls="dashed", lw=3, label="Fault plane 2")
        # ax1[i].axvline(
        #     x=mean_Full2[i] + std_Full2[i], c="steelblue", lw=1, ls="--", label="std 2", alpha=0.5
        # )
        # ax1[i].axvline(x=mean_Full2[i] - std_Full2[i], c="steelblue", lw=1, ls="--", alpha=0.5)
        # ax1[i].axvline(x=f_planes[i], c="steelblue", lw=1, ls="--", alpha=0.5)
        # ax1[i].axvline(x=aux_planes[i], c="steelblue", lw=1, ls="--", alpha=0.5)

        ax1[i].set_xlabel(sdr_names[i], fontsize=45)
        if i == 0:
            ax1[i].set_ylabel("Probability", fontsize=45)

        ax1[i].tick_params(axis="x", labelsize=23)
        ax1[i].tick_params(axis="y", labelsize=23)
        ax1[i].ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
        tx = ax1[i].yaxis.get_offset_text()
        tx.set_fontsize(23)
        ax1[i].yaxis.label.set_color("blue")
        ax1[i].tick_params(axis="y", colors="blue")
        ax1[i].xaxis.set_ticks(ticklabels[i])
        ax1[i].set_xticklabels(ticklabels[i])
        ax1[i].locator_params(axis="y", nbins=4)

        if event_name == "S0183a":
            pass
        else:
            if i == 2:
                ax2.set_ylabel("Probability", fontsize=45)
            ax2.locator_params(axis="y", nbins=4)
            ax2.yaxis.label.set_color("red")
            ax2.tick_params(axis="y", colors="red")
            ax2.ticklabel_format(style="sci", axis="y", scilimits=(-2, 2))
            tx = ax2.yaxis.get_offset_text()
            tx.set_fontsize(23)
            ax2.tick_params(axis="y", labelsize=23)
            if i == 2:
                # where some data has already been plotted to ax
                handles, labels = ax1[i].get_legend_handles_labels()

                # manually define a new patch
                patch = mpatches.Patch(color="red", label="Direct", alpha=0.8)

                # handles is a list, so append manual patch
                handles.append(patch)

                # ax1[i].legend(handles=handles, ncol=2)
                ax1[i].legend(handles=handles, ncol=2)

        ax1[i].set_xlim(sdr_mins[i], sdr_maxs[i])
        # ax1[i].set_ylim(0, 0.01)

    fig1.suptitle(event_name, fontsize=30)

    fig = Plot_GS_BB(
        MT_sdrs[:, 0],
        MT_sdrs[:, 1],
        MT_sdrs[:, 2],
        azimuths=[10, 10, 10],
        inc_angles=[10, 10, 10],
        phase_names=["P", "P", "P"],
        color="blue",
    )
    plt.savefig(
        pjoin(h5_file_folder, f"Overal_BBB__{event_name}.svg",), dpi=300,
    )
    plt.close()
    return fig1


def polarization(
    h5_file_folder: str,
    event: obspy.core.event.Event,
    method: str,
    misfit_name: str,
    fwd: _Forward._AbstractForward,
    phases: [str],
    components: [str],
    depths: [float],
    DOF: float,
    fmin: float = None,
    fmax: float = None,
):
    amount_of_phases = len(phases)
    n_lowest = 1
    for idepth, depth in enumerate(depths):
        print(depth)
        GS_File = glob.glob(
            pjoin(
                h5_file_folder,
                f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )
        )[0]
        Direct_File = glob.glob(
            pjoin(
                h5_file_folder,
                f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.hdf5",
            )
        )[0]

        ## ================ READ GS =============================
        (depth_GS, sdr, M0_GS, misfit_L2_GS,) = _ReadH5.Read_GS_h5(
            Filename=GS_File, amount_of_phases=amount_of_phases
        )
        Total_L2_GS = np.sum(misfit_L2_GS, axis=1)
        # Total_L2_norm_GS = np.sum(misfit_L2_norm_GS, axis=1)
        # GOF = ( (Total_L2_GS - DOF ) * 100 ) / DOF
        # GOF_GS = (Total_L2_norm_GS / DOF) * 100
        GOF_GS = Total_L2_GS / DOF

        lowest_indices = Total_L2_GS.argsort()[0:n_lowest]
        # lowest_indices = GOF_GS.argsort()[0:n_lowest]

        sdr = sdr[lowest_indices, :]
        print("strike", sdr[0][0], "dip", sdr[0][1], "rake", sdr[0][2])
        depth_GS = depth_GS[lowest_indices]
        M0_GS = M0_GS[lowest_indices]
        # shifts["P"] = shifts["P"][lowest_indices]
        # shifts["S"] = shifts["S"][lowest_indices]

        # L2_GS = np.append(L2_GS, Total_L2_GS[lowest_indices][0])
        # L2_GS = np.append(L2_GS, GOF_GS[lowest_indices][0])
        # if depth == 8:
        #     lowest_indices[0] = lowest_indices[2]
        #     L2_GS = np.append(L2_GS, GOF_GS[lowest_indices][0])
        # else:
        #     L2_GS = np.append(L2_GS, GOF_GS[lowest_indices][0])

        ## ============ Read Direct ========================
        (
            depth_Direct,
            FULL_MT,
            DC_MT,
            CLVD_MT,
            misfit_L2_Direct,
            Epsilon,
            M0,
            M0_DC,
            M0_CLVD,
            angles,
            cond_nr,
        ) = _ReadH5.Read_Direct_Inversion(Direct_File, amount_of_phases=amount_of_phases)

        Total_L2_Direct = np.sum(misfit_L2_Direct)
        # Total_L2_norm_Direct = np.sum(misfit_L2_norm_Direct, axis=1)
        # GOF_Direct = (Total_L2_norm_Direct / DOF) * 100
        GOF_Direct = Total_L2_Direct / DOF
        # L2_Direct = np.append(L2_Direct, Total_L2_Direct[0])
        # L2_Direct = np.append(L2_Direct, GOF_Direct)

        RP = np.array([])
        RSV = np.array([])
        RSH = np.array([])
        incs = np.array([])
        take_off_with_az = {}
        az = event.az

        azs = np.arange(0, 360, 10)
        to_angles = np.arange(0, 90, 5)

        Combs = np.array(np.meshgrid(azs, to_angles)).T.reshape(-1, 2)
        azs = Combs[:, 0]
        # distances = Combs[:,1]
        angles = Combs[:, 1]

        strike = sdr[0][0]
        dip = sdr[0][1]
        rake = sdr[0][2]
        phase_names = ["P", "S"]
        for j, comb in enumerate(Combs):
            az = comb[0]
            to_angle = comb[1]

            inc = to_angle
            incs = np.append(incs, inc)
            RP = np.append(
                RP,
                _RadiationPattern.R_P(
                    take_off_angle=inc, strike=strike, dip=dip, rake=rake, az=az
                ),
            )
            RSV = np.append(
                RSV,
                _RadiationPattern.R_SV(
                    take_off_angle=inc, strike=strike, dip=dip, rake=rake, az=az
                ),
            )
            RSH = np.append(
                RSH,
                _RadiationPattern.R_SH(
                    take_off_angle=inc, strike=strike, dip=dip, rake=rake, az=az
                ),
            )
        RP.shape = (RP.size, 1)
        RSV.shape = (RSV.size, 1)
        RSH.shape = (RSH.size, 1)

        tt = fwd.taup_veloc.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=event.distance, phase_list=["P", "S"]
        )
        take_off_with_az = {}
        take_off_with_az[tt[0].name] = [tt[0].takeoff_angle, event.az]
        take_off_with_az[tt[1].name] = [tt[1].takeoff_angle, event.az]

        fig = Plot_Radiation_BB(
            [strike, dip, rake],
            RP,
            RSV,
            RSH,
            azs,
            incs,
            "grey",
            take_offs_with_azs=take_off_with_az,
            plot_take_offs=True,
        )
        plt.savefig(
            pjoin(
                h5_file_folder,
                f"Polarization_{event.name}_{depth}_{fmin}_{fmax}_{misfit_name}_{fwd.veloc_name}.svg",
            ),
            dpi=300,
        )
        plt.close()


def Plot_Radiation_BB(
    MT,
    RP,
    RSV,
    RSH,
    azimuths,
    inc_angles,
    color,
    height=None,
    take_offs_with_azs=None,
    plot_take_offs=True,
    horizontal=False,
):

    if horizontal:
        width = 15.0
        height = 7.0

        axis_heigt = 5.0 / height
        resid_heigt = 1.0 - axis_heigt
        title_height = resid_heigt

        axis_width = 5.0 / width
    else:
        if height == None:
            height = 19.0
        axis_height = 5.0 / height
        resid_height = 1.0 - 3.0 * axis_height
        title_height = resid_height / 3.0

    ## Full moment tensor:
    img1, buf1 = Get_bb_img(MT, color, alpha=0.4)

    if horizontal:
        fig = plt.figure(figsize=(width, height), dpi=200)
        ax_1 = fig.add_axes([2 * (axis_width), 0.0, axis_width, axis_heigt])
    else:
        fig = plt.figure(figsize=(5, height), dpi=200)
        ax_1 = fig.add_axes([0.0, 2 * (axis_height + title_height), 1.0, axis_height])

    if img1 is not None:
        ax_1.imshow(img1 / np.max(img1.flatten()))
    if horizontal:
        ax_X = fig.add_axes([2 * (axis_width), 0.0, axis_width, axis_heigt], label="Circle_ray")
    else:
        ax_X = fig.add_axes(
            [0.0, 2 * (axis_height + title_height), 1.0, axis_height], label="Circle_ray"
        )
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)

    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(RP)
    x_scaled = 0.5 * (RP) + 0.5  # min_max_scaler.fit_transform(RSV)
    df_S = pd.DataFrame(x_scaled)
    c2 = df_S[0]
    colors = [cm.seismic(color) for color in c2]

    x_tot = np.array([])
    y_tot = np.array([])

    if azimuths is not None and inc_angles is not None:
        for a, i, color_i in zip(azimuths, inc_angles, colors):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0

            x_tot = np.append(x_tot, x)
            y_tot = np.append(y_tot, y)

            p = Circle(
                (x, y),
                0.02,
                linewidth=2,
                edgecolor=color_i,
                color=color_i,
                zorder=0,
                facecolor="k",
                fill=True,
            )
            ax_X.add_patch(p)

    if take_offs_with_azs is not None:
        for phase in take_offs_with_azs:
            a = take_offs_with_azs[phase][1]
            i = take_offs_with_azs[phase][0]
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            if plot_take_offs:
                p = Circle(
                    (x, y), 0.02, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
                )
                ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=25)
    for a in [ax_1, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")

    sc = plt.scatter(
        x_tot, y_tot, s=0, c=RP[:, 0], vmin=-1.0, vmax=1.0, cmap="seismic", facecolors="none",
    )
    if horizontal:
        cbar_1 = fig.add_axes([2 * axis_width, axis_heigt, axis_width, title_height])
    else:
        cbar_1 = fig.add_axes([0.0, 3 * axis_height + 2 * title_height, 1.0, title_height])
    cbar_1.set_xticks([])
    cbar_1.set_yticks([])
    cbar_1.axis("off")
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cbaxes = inset_axes(cbar_1, width="80%", height="25%", loc="center")
    cbar = plt.colorbar(sc, cax=cbaxes, orientation="horizontal")
    cbar.set_label("Radiation P", rotation=0, labelpad=10, fontsize=25)

    ###
    if horizontal:
        ax_2 = fig.add_axes([axis_width, 0.0, axis_width, axis_heigt])
    else:
        ax_2 = fig.add_axes([0.0, axis_height + title_height, 1.0, axis_height])
    img2, buf2 = Get_bb_img(MT, color, alpha=0.4)

    if img2 is not None:
        ax_2.imshow(img2 / np.max(img2.flatten()))

    if horizontal:
        ax_X = fig.add_axes([axis_width, 0.0, axis_width, axis_heigt], label="Circle_ray")
    else:
        ax_X = fig.add_axes(
            [0.0, axis_height + title_height, 1.0, axis_height], label="Circle_ray"
        )
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = 0.5 * (RSV) + 0.5  # min_max_scaler.fit_transform(RSV)
    df_S = pd.DataFrame(x_scaled)
    c2 = df_S[0]
    colors = [cm.seismic(color) for color in c2]

    x_tot = np.array([])
    y_tot = np.array([])

    if azimuths is not None and inc_angles is not None:
        for a, i, color_i in zip(azimuths, inc_angles, colors):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0

            x_tot = np.append(x_tot, x)
            y_tot = np.append(y_tot, y)

            p = Circle(
                (x, y),
                0.02,
                linewidth=2,
                edgecolor=color_i,
                color=color_i,
                zorder=0,
                facecolor="k",
                fill=True,
            )
            ax_X.add_patch(p)
    if take_offs_with_azs is not None:
        for phase in take_offs_with_azs:
            a = take_offs_with_azs[phase][1]
            i = take_offs_with_azs[phase][0]
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            if plot_take_offs:
                p = Circle(
                    (x, y), 0.02, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
                )
                ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=25)

    for a in [ax_2, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")

    sc = plt.scatter(
        x_tot,
        y_tot,
        s=0,
        c=RSV[:, 0],
        vmin=-1.0,  # RP.min(),
        vmax=1.0,  # RP.max(),
        cmap="seismic",
        facecolors="none",
    )

    if horizontal:
        cbar_2 = fig.add_axes([axis_width, axis_heigt, axis_width, title_height])
    else:
        cbar_2 = fig.add_axes([0.0, 2 * axis_height + title_height, 1.0, title_height])
    cbar_2.set_xticks([])
    cbar_2.set_yticks([])
    cbar_2.axis("off")
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cbaxes = inset_axes(cbar_2, width="80%", height="25%", loc="center")
    cbar = plt.colorbar(sc, cax=cbaxes, orientation="horizontal")
    cbar.set_label("Radiation SV", rotation=0, labelpad=10, fontsize=25)

    ### Radiation SH

    if horizontal:
        ax_3 = fig.add_axes([0.0, 0.0, axis_width, axis_heigt])
    else:
        ax_3 = fig.add_axes([0.0, 0.0, 1.0, axis_height])
    img3, buf3 = Get_bb_img(MT, color, alpha=0.4)

    if img3 is not None:
        ax_3.imshow(img3 / np.max(img3.flatten()))

    if horizontal:
        ax_X = fig.add_axes([0.0, 0.0, axis_width, axis_heigt], label="Circle_ray")
    else:
        ax_X = fig.add_axes([0.0, 0.0, 1.0, axis_height], label="Circle_ray")
    ax_X.set_xlim((-1, 1))
    ax_X.set_ylim((-1, 1))
    p = Circle((0.0, 0,), 0.99, linewidth=2, edgecolor="k", zorder=0, fill=False)
    ax_X.add_patch(p)
    # min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = 0.5 * (RSH) + 0.5  # min_max_scaler.fit_transform(RSV)
    # x_scaled = min_max_scaler.fit_transform(RSH)
    df_S = pd.DataFrame(x_scaled)
    c2 = df_S[0]
    colors = [cm.seismic(color) for color in c2]

    x_tot = np.array([])
    y_tot = np.array([])

    if azimuths is not None and inc_angles is not None:
        for a, i, color_i in zip(azimuths, inc_angles, colors):
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0

            x_tot = np.append(x_tot, x)
            y_tot = np.append(y_tot, y)

            p = Circle(
                (x, y),
                0.02,
                linewidth=2,
                edgecolor=color_i,
                color=color_i,
                zorder=0,
                facecolor="k",
                fill=True,
            )
            ax_X.add_patch(p)
    if take_offs_with_azs is not None:
        for iphase, phase in enumerate(take_offs_with_azs):
            a = take_offs_with_azs[phase][1]
            i = take_offs_with_azs[phase][0]
            if i > 90.0:
                x = np.sin(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
                y = np.cos(np.deg2rad(a + 180)) * (180.0 - i) / 90.0
            else:
                x = np.sin(np.deg2rad(a)) * i / 90.0
                y = np.cos(np.deg2rad(a)) * i / 90.0
            if plot_take_offs:
                p = Circle(
                    (x, y), 0.02, linewidth=2, edgecolor="k", zorder=0, facecolor="k", fill=True
                )
                ax_X.add_patch(p)
            ax_X.text(x - 0.005, y + 0.03, s=phase, fontsize=25)

    for a in [ax_3, ax_X]:
        a.set_xticks([])
        a.set_yticks([])
        a.axis("off")

    sc = plt.scatter(
        x_tot, y_tot, s=0, c=RSH[:, 0], vmin=-1.0, vmax=1.0, cmap="seismic", facecolors="none",
    )
    if horizontal:
        cbar_3 = fig.add_axes([0.0, axis_heigt, axis_width, title_height])
    else:
        cbar_3 = fig.add_axes([0.0, axis_height, 1.0, title_height])
    cbar_3.set_xticks([])
    cbar_3.set_yticks([])
    cbar_3.axis("off")
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cbaxes = inset_axes(cbar_3, width="80%", height="25%", loc="center")
    cbar = plt.colorbar(sc, cax=cbaxes, orientation="horizontal")
    cbar.set_label("Radiation SH", rotation=0, labelpad=10, fontsize=25)

    return fig


def k_mean_distance(data, cx, i_centroid, cluster_labels):
    distances = [np.sqrt((x - cx) ** 2) for x in data[cluster_labels == i_centroid]]
    return np.mean(distances)


def k_mean_distance_2d(data, cx, cy, i_centroid, cluster_labels):
    distances = [
        np.sqrt((x - cx) ** 2 + (y - cy) ** 2) for (x, y) in data[cluster_labels == i_centroid]
    ]
    return np.mean(distances)


def k_mean_distance_3d(data, cx, cy, cz, i_centroid, cluster_labels):
    distances = [
        np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        for (x, y, z) in data[cluster_labels == i_centroid]
    ]
    return np.mean(distances)
