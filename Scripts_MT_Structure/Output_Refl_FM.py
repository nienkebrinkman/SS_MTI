from os import listdir
from os.path import isfile, join
import obspy
import matplotlib.pyplot as plt
import numpy as np


folder = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_BKE/mxx_new/"

st_files = [f for f in listdir(folder) if f.startswith("st") if isfile(join(folder, f))]

st = obspy.Stream()
dists = []
azis = []
for st_file in st_files:
    st_temp = obspy.read(join(folder, st_file))
    st_temp[0].stats.channel = "XB" + st_file.split(".")[-1]
    distance = st_temp[0].stats.sac.gcarc
    st_temp[0].stats.distance = distance

    B = st_temp[0].stats.sac.b  # beginning time
    tstart = st_temp[0].stats.starttime  # absolute starttime
    st_temp[0] = st_temp[0].trim(
        starttime=tstart + B, endtime=st_temp[0].stats.endtime, pad=True, fill_value=0.0
    )
    # st_temp[0].stats.starttime = tstart + B  # absolute time of event

    st += st_temp
    dists.append(distance)
    azis.append(st_temp[0].stats.sac.az)

# st_select = st.select(channel="XB" + "Z")
# st_select.plot(type="section")


def filter_tr(tr, fmin=1.0 / 10.0, fmax=1.0 / 2, zerophase=False):
    tr.filter("highpass", freq=fmin, corners=4, zerophase=zerophase)
    tr.filter("highpass", freq=fmin, corners=4, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, corners=4, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, corners=4, zerophase=zerophase)


st[0].stats.channel = st[0].stats.channel.replace("t", "T")
st[1].stats.channel = st[1].stats.channel.replace("r", "R")
st[2].stats.channel = st[2].stats.channel.replace("z", "Z")

# for tr in st:
#     filter_tr(tr, fmin=0.1, fmax=0.9, zerophase=False)
st = st.rotate(method="RT->NE", back_azimuth=90.0)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(28, 12))

colours = ["r", "b", "k"]

for i, comp in enumerate(["Z", "N", "E"]):
    st_select = st.select(channel="XB" + comp)
    print(st_select[0].stats.starttime)
    for j, tr in enumerate(st_select):
        P_arr = ((89.3 / np.cos(np.deg2rad(dists[j]))) / 4.8036) / tr.stats.delta
        ax.plot(
            tr.times(),
            tr.data + j * 1e-9,
            colours[i],
            label=f"Comp: {comp}, dist: {dists[j]}, azi: {azis[j]} ",
        )
        # ax.plot([P_arr, P_arr], [-0.5e-9 + j * 1e-9, 0.5e-9 + j * 1e-9])
plt.legend()
plt.show()

