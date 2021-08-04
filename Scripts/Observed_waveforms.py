import obspy
import numpy as np
import matplotlib.pyplot as plt
from os.path import join as pjoin

from argparse import ArgumentParser
from mqs_reports.catalog import Catalog
from obspy import UTCDateTime as utct

import SS_MTI


def filter(tr, fmin=1.0 / 10.0, fmax=1.0 / 2, zerophase=True):
    tr.filter("highpass", freq=fmin, zerophase=zerophase)
    tr.filter("highpass", freq=fmin, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, zerophase=zerophase)


save_path = "/home/nienke/Documents/Research/Data/MTI/"

path = "/home/nienke/Documents/Research/Data/MTI/Catalog/"
path_to_inventory = pjoin(path, "inventory.xml")
path_to_catalog = pjoin(path, "catalog.xml")
inv = SS_MTI.DataGetter.read_inv(inv_path=path_to_inventory)  # Inventory file
cat = SS_MTI.DataGetter.read_cat(cat_path=path_to_catalog)  # Catalog file

phases = 2
# fig, ax = plt.subplots(nrows=phases, ncols=1, sharex="all", figsize=(28, 8))
fig, ax = plt.subplots(nrows=phases, ncols=1, sharex="all", figsize=(25, 18))
# Fig 3a
Ylims = {"P": [0.8e-9, 0.8e-9], "S": [2.5e-9, 2.5e-9, 2.5e-9]}
# Fig 3b
# Ylims = {'P': [1.2e-9, 1.2e-9], 'S': [3.8e-9, 3.8e-9, 3.8e-9]}

obs_data = obspy.Stream()
obs_data_long = obspy.Stream()

offset_plot = False

for i, source in enumerate(np.array([0, 1, 2, 3, 4])):
    if source == 0:
        ## S0235b
        event = cat.select(name="S0235b").events[0]
        event.read_waveforms(inv=inv, sc3dir="/mnt/marshost/sc3data/")
        phases = ["P", "S", "S", "P", "S"]
        components = ["Z", "T", "Z", "R", "R"]
        corrs = [0.2 + 1.24, 10.4 - 9.2, 11.1 - 9.2, 0.2 + 1.24, 11.1 - 9.2]
        # # phases = ['P', 'P']
        # # components = ['Z', 'R']
        # # corrs = [0., 0.]
        #
        # phases = ['S', 'S', 'S']
        # components = ['Z', 'R', 'T']
        # corrs = [10.1, 10.1,10.1]
        # Ylims = {'P': [2e-9, 2e-9], 'S': [4e-9, 4e-9, 4e-9]}  # ZRT
        color = "red"
        fmin = 0.1  # 1. / 10.#1./8.
        fmax = 0.5  # 1. / 1.#1./5.

    elif source == 1:
        ## S0173a
        event = cat.select(name="S0173a").events[0]
        event.read_waveforms(inv=inv, sc3dir="/mnt/marshost/sc3data/")
        phases = ["P", "S", "S", "P", "S"]
        components = ["Z", "T", "Z", "R", "R"]
        # corrs = [-0.8,2.,2.,-0.8,2.]
        corrs = [-0.3, 2.9, 2.0, -0.3, 2.9]
        # phases = ['P', 'P']
        # components = ['Z', 'R']
        # corrs = [-0.8, -0.8]

        # phases = ['S', 'S', 'S']
        # components = ['Z', 'R', 'T']
        # corrs = [2., 2.,2.]
        # Ylims = {'P': [2e-9, 2e-9], 'S': [4e-9, 4e-9, 4e-9]}  # ZRT
        color = "green"
        fmin = 0.1  # 1. / 10.#1./8.
        fmax = 0.4  # 1. / 2.#1./5.

    elif source == 2:
        ## S0183a
        event = cat.select(name="S0183a").events[0]
        event.read_waveforms(inv=inv, sc3dir="/mnt/marshost/sc3data/")
        event.distance = 44.5
        event.baz = 73
        event.az = 253
        event.latitude = 15.09
        event.longitude = 179.59

        phases = ["P", "S", "S", "P", "S"]
        components = ["Z", "T", "Z", "R", "R"]
        # corrs = [2.1,0.,0.,2.1,0.]
        corrs = [2.1 - 3.17, 2.2 - 7.075, 1.2 - 7.075, 2.1 - 3.17, 1.2 - 7.075]

        # phases = ['P', 'P']
        # components = ['Z', 'R']
        # corrs = [2.1, 2.1]
        # Ylims = {'P': [4e-10, 4e-10], 'S': [5e-10, 1e-9, 2e-9]}  # ZRT
        color = "blue"
        fmin = 1.0 / 5.0  # 1. / 10.#1./8.
        fmax = 1.0 / 2.5  # 1. / 2.#1./5.
    elif source == 3:
        ## S0809a:
        event = cat.select(name="S0809a").events[0]
        event.read_waveforms(inv=inv, sc3dir="/mnt/marshost/sc3data/")
        phases = ["P", "S", "S", "P", "S"]
        components = ["Z", "T", "Z", "R", "R"]
        corrs = [1.2, 0.2, 0.2, 1.2, 0.2]
        # corrs = [-0.7, 0.2, 0.9, -0.7, 0.9]
        # corrs = [0.0, 0.0, 0.0, 0.0, 0.0]
        # tstars= [0.4, 0.2, 0.2, 0.4, 0.2]
        fmin = 0.2
        fmax = 0.8
        color = "purple"
    elif source == 4:
        ## S0820a:
        event = cat.select(name="S0820a").events[0]
        event.read_waveforms(inv=inv, sc3dir="/mnt/marshost/sc3data/")
        phases = ["P", "S", "S", "P", "S"]
        components = ["Z", "T", "Z", "R", "R"]
        corrs = [8.5, 2.8, 2.8, 8.5, 2.8]
        fmin = 0.2
        fmax = 0.8
        color = "black"

    print(f"{event.name}, P:")
    print(event.picks["P"])
    print(f"{event.name}, S:")
    print(event.picks["S"])
    print(f"lat: {event.latitude}")
    print(f"lon: {event.longitude}")
    # For plot 3a:
    # fmin = 1.0 / 5.0  # 1. / 10.#1./8.
    # fmax = 1.0 / 2.5  # 1. / 2.#1./5.

    fmin = 0.2
    fmax = 0.4

    nphase = 2  # len(phases)
    arrs_real = []
    for iphase in range(0, nphase):
        arrs_real.append(utct(event.picks[phases[iphase]]) - event.origin_time + corrs[iphase])

    normfac = np.zeros(nphase)
    st_real = obspy.Stream()
    st_real_orig = obspy.Stream()
    st_sigma = obspy.Stream()

    normfac_phase = {"P": None, "S": None}
    LQT_value = False
    for iphase in range(nphase):
        # if source == 2:
        #     if iphase == 1 or iphase == 2 or iphase == 4:
        #         continue
        if components[iphase] in ["Z", "N", "E"]:
            tr_orig = event.waveforms_VBB.select(channel="BH" + components[iphase]).copy()[0]
        else:
            st_raw = event.waveforms_VBB.copy()
            st_raw.rotate(method="NE->RT", back_azimuth=event.baz)
            tr_orig = st_raw.select(channel="BH" + components[iphase])[0]

        filter(tr_orig, fmin=fmin, fmax=fmax, zerophase=False)
        arrs = arrs_real
        if source == 1 and phases[iphase] == "P":
            tr_real = tr_orig.slice(
                starttime=event.origin_time + arrs[iphase] - 1.0,
                endtime=event.origin_time + arrs[iphase] + 17.0,
            )
            ax[iphase].axvspan(17, 40, facecolor="green", alpha=0.1)
            y = ax[iphase].get_ylim()[0] * 0.8
            ax[iphase].text(
                35, y, "Glitch", verticalalignment="center", color="green", fontsize=35
            )
            ax[iphase].axvline(x=17, c="green", ls="dotted", alpha=0.7)
            ax[iphase].axvline(x=40, c="green", ls="dotted", alpha=0.7)
        else:
            tr_real = tr_orig.slice(
                starttime=event.origin_time + arrs[iphase] - 1.0,
                endtime=event.origin_time + arrs[iphase] + 30.0,
            )
        tr_orig.trim(
            starttime=event.origin_time + arrs[iphase] - 60.0,
            endtime=event.origin_time + arrs[iphase] + 120.0,
        )
        tr_sigma = tr_orig.slice(
            starttime=event.origin_time + arrs[iphase] - 30.0,
            endtime=event.origin_time + arrs[iphase] - 20.0,
        )
        Sigma = np.std(tr_sigma.data)

        # SCALING THE DATA
        if source == 0:
            obs_data += tr_real
            obs_data_long += tr_orig
        else:
            dt = tr_orig.stats.delta
            # For plot 3b: comment the multiplication with M0_corr
            M0_corr = (np.max(np.abs(obs_data[iphase].data[: int(10.0 / dt)]))) / (
                np.max(np.abs(tr_real.data[: int(10.0 / dt)]))
            )
            tr_real.data = tr_real.data * M0_corr
            tr_orig.data = tr_orig.data * M0_corr

        normfac[iphase] = max(abs(tr_orig.data))
        if (normfac_phase[phases[iphase]] is None) and (components[iphase] == "Z"):
            normfac_phase[phases[iphase]] = normfac[iphase]
        elif (normfac_phase[phases[iphase]] is None) and (components[iphase] == "T"):
            normfac_phase[phases[iphase]] = normfac[iphase]
        if offset_plot:
            if phases[iphase] == "P":
                limit = Ylims["P"][0]
            elif phases[iphase] == "S":
                limit = Ylims["S"][0] * 1.5
            else:
                limit = 0.0
        else:
            limit = 0.0
        ax[iphase].plot(
            tr_real.times() - 1.0, tr_real.data + (i * limit), color=color, lw=4, label=event.name,
        )
        ax[iphase].plot(tr_orig.times() - 60.0, tr_orig.data + (i * limit), color=color, lw=0.5)
        if source == 0:
            ax[iphase].axvspan(-1, 10, facecolor="grey", alpha=0.2)
        # if source ==1:
        #     Weight_vector = np.linspace(1,3, len(tr_real.data))
        #     Weight_vector[:int(7.0 / 0.05)] = 1
        #     Weight_vector[int(7.0 / 0.05):] = 3
        #     ax[iphase].plot(tr_real.times() - 1.,
        #                     Sigma * Weight_vector,
        #                     lw=0.8, ls='dotted', c='grey',
        #                     label=r'$Weight$''\n'r'$(\sigma~based)$')

        if i == 0:
            ax[iphase].text(
                s="%s%s" % (phases[iphase], components[iphase]),
                x=0.98,
                y=0.9,  # 0.72,
                ha="right",
                transform=ax[iphase].transAxes,
                color="black",
                fontsize=45,
            )

        ax[iphase].tick_params(axis="both", which="major", labelsize=35)
        ax[iphase].tick_params(axis="both", which="minor", labelsize=25)
        if not offset_plot:
            if phases[iphase] == "P":
                if components[iphase] == "Z" or components[iphase] == "L":
                    ax[iphase].set_ylim(-Ylims["P"][0], Ylims["P"][0])
                elif components[iphase] == "R" or components[iphase] == "Q":
                    ax[iphase].set_ylim(-Ylims["P"][1], Ylims["P"][1])
                # Fig3a:
                ax[iphase].yaxis.set_ticks([-Ylims["P"][0], 0, Ylims["P"][0]])
                ax[iphase].set_yticklabels([-1, 0, 1])

            else:
                if components[iphase] == "Z" or components[iphase] == "L":
                    ax[iphase].set_ylim(-Ylims["S"][0], Ylims["S"][0])
                elif components[iphase] == "R" or components[iphase] == "Q":
                    ax[iphase].set_ylim(-Ylims["S"][1], Ylims["S"][1])
                elif components[iphase] == "T":
                    ax[iphase].set_ylim(-Ylims["S"][2], Ylims["S"][2])
                # Fig3a:
                ax[iphase].yaxis.set_ticks([-Ylims["S"][2], 0, Ylims["S"][2]])
                ax[iphase].set_yticklabels([-1, 0, 1])
        else:
            ax[iphase].yaxis.set_ticklabels([])
# ax[-1].set_xlim(-10, 60.0)
ax[-1].set_xlim(-10, 31.0)
# fig.text(0.45, 0.88, "Observed waveforms", ha="center",
#          va="bottom", size="x-large", color="blue")
fig.text(0.01, 0.5, "Normalized \n Displacement", va="center", rotation="vertical", fontsize=45)
# fig.text(0.01, 0.5, 'Displacement (nm)',
#          va='center', rotation='vertical', fontsize=45)
# ax[iphase].ticklabel_format(
#     style="sci", axis='y', scilimits=(-2, 2))

for iphase in range(nphase):
    ax[iphase].axvline(x=0.0, c="black")
    # ax[iphase].axvline(x=30., c='black')
    ax[iphase].axvline(x=-1.0, c="grey", ls="dashed")
    if phases[iphase] == "P" and source == 1:
        ax[iphase].axvline(x=17.0, c="grey", ls="dashed")
    else:
        ax[iphase].axvline(x=30.0, c="grey", ls="dashed")
    ax[iphase].axvline(x=0.0, c="black")

    if phases[iphase] == "P":
        if components[iphase] == "Z":
            ax[0].legend(prop={"size": 14}, loc="center left", bbox_to_anchor=(0.0, 1.1))
        y = ax[iphase].get_ylim()[1] * 0.85
        ax[iphase].text(0.1, y, "P", verticalalignment="center", color="black", fontsize=35)
    elif phases[iphase] == "S":
        y = ax[iphase].get_ylim()[1] * 0.85
        ax[iphase].text(
            0.1, y, phases[iphase], verticalalignment="center", color="black", fontsize=35
        )
    ax[iphase].get_yaxis().get_offset_text().set_visible(False)
    ax_max = max(ax[iphase].get_yticks())
    exponent_axis = np.floor(np.log10(ax_max)).astype(int)
    # ax[iphase].annotate(r'$\times$10$^{%i}$' % (exponent_axis),
    #                     xy=(.01, .82), xycoords='axes fraction', fontsize = 35)

ax[-1].set_xlabel("time after phase (s)", fontsize=45)
# plt.tight_layout()
# plt.show()
plt.savefig(save_path + "Observed_Waveforms.pdf", dpi=300)
plt.close()

from mpl_toolkits.basemap import Basemap

mars_dir = "/home/nienke/Documents/Research/Data/mars_pictures/Mars_lightgray.jpg"

fig = plt.figure(figsize=(10, 8))

# m = Basemap(projection='moll', lon_0=round(0.0))
m = Basemap(
    projection="merc", llcrnrlat=-40, urcrnrlat=40, llcrnrlon=120, urcrnrlon=200, resolution="c",
)

# draw parallels and meridians.
par = np.arange(-90, 90, 30)
label_par = np.full(len(par), True, dtype=bool)
meridians = np.arange(-180, 180, 30)
label_meri = np.full(len(meridians), True, dtype=bool)

m.drawmeridians(np.arange(-180, 180, 30), labels=label_meri)
m.drawparallels(np.arange(-90, 90, 30), label=label_par)

m.warpimage(mars_dir)
lat_rec = 4.5  # 02384
lon_rec = 135.623447
mstatlon, mstatlat = m(lon_rec, lat_rec)
m.plot(mstatlon, mstatlat, "k^", markersize=20, label="InSight")

EQlonA, EQlatA = m(162.8890468, 11.40707267)  # S0235b
EQlonB, EQlatB = m(164.9642605, 3.437219264)  # S0173a
EQlonC, EQlatC = m(179.59, 15.09)  # S0183a
EQlonD, EQlatD = m(164.6489946, 5.390502842)  # S0809a
EQlonE, EQlatE = m(165.3877124, -1.629058007)  # S0820a
m.plot(EQlonA, EQlatA, "r*", markersize=20, zorder=10, label="S0235b")
m.plot(EQlonB, EQlatB, "g*", markersize=20, zorder=10, label="S0173a")
m.plot(EQlonC, EQlatC, "b*", markersize=20, zorder=10, label="S0183a")
m.plot(
    EQlonD,
    EQlatD,
    color="purple",
    marker="*",
    linestyle="None",
    markersize=20,
    zorder=10,
    label="S0809a",
)
m.plot(EQlonE, EQlatE, "k*", markersize=20, zorder=10, label="S0820a")
plt.legend(fontsize=20)
plt.tight_layout()
plt.savefig(save_path + "location_events.png", dpi=300)
plt.close()
