""" Converting output file to hdf5 """
import h5py
from os import listdir
from os.path import isfile, join
import numpy as np
import obspy


# folder = "/home/nienke/Documents/Research/SS_MTI/External_packages/erzsol3/ERZSOL3/"
folder = "/home/nienke/Documents/Research/SS_MTI/External_packages/erzsol3/Test_Files/"

or_time = obspy.UTCDateTime(2019, 1, 1, 12, 0, 0)

erz_folder = folder
erzFiles = [f for f in listdir(erz_folder) if f.endswith(".tx.z") if isfile(join(erz_folder, f))]

filenameH5Output = join(folder, "Test.h5")

n_examples = len(erzFiles)
n_clusters = 1
# ns = 4096s
cmd_folder = folder

sou_physics = np.zeros((n_examples, 4), dtype="int")  # strike, dip, rake, seismic moment
sou_coordinates = np.zeros((n_examples, 3), dtype="int")
one_hot_vectors = np.zeros((n_examples, n_clusters), dtype="int")
azimuths = np.zeros((n_examples, 1))
ranges = np.zeros((n_examples, 1))


# ncomp = 0
# # Read single file to get number of receivers to initialize data_matrix
# f = open(join(erz_folder, erzFiles[0]), "rb")
# k = 4
# f.seek(k)
# nt = np.fromfile(f, dtype="int32", count=1)[0]  # number of receivers
# data_matrix = np.zeros((n_examples, nt, ns))
# f.close()

# # Begin loop over all the input files
# for i, ef in enumerate(erzFiles):

#     f = open(join(erz_folder, ef), "rb")

#     cmd_file = join(cmd_folder, ef.split(".")[0] + ".cmd")
#     f_cmd = open(cmd_file, "r")

#     lines = f_cmd.read().splitlines()
#     f_cmd.close()

#     # sou_physics[i, :] = np.fromstring(lines[28], dtype="int", sep=" ")
#     # sou_coordinates[i, :] = np.fromstring(lines[31], dtype="int", sep=" ")
#     # one_hot_vectors[i, :] = np.fromstring(lines[34], dtype="int", sep=" ")

#     # First part information about number of receivers and components per receiver
#     k = 4
#     f.seek(k)
#     n_rec = np.fromfile(f, dtype="int32", count=1)[0]  # number of receivers
#     k += 4
#     f.seek(k)
#     n_comp = np.fromfile(f, dtype="int32", count=1)[0]  # number of components per receiver
#     k += 8

#     # Not best prgramming. But does the job. These are the bytes at which to read beta info and data:
#     num_bytes = np.array([4, 4, 5, 3, 4, 4, 4, 4, 4, ns * 4, 4])

#     # Loop over all the receivers and their individual components
#     for i_r in range(0, n_rec):
#         for j in range(0, n_comp):

#             k += num_bytes[0]
#             f.seek(k)
#             dist = np.fromfile(f, dtype="float32", count=1)
#             k += num_bytes[1]
#             f.seek(k)
#             azi = np.fromfile(f, dtype="float32", count=1)
#             k += num_bytes[2]
#             f.seek(k)
#             comp = np.fromfile(f, dtype="|S1", count=1).astype(str)[0]
#             k += num_bytes[3]
#             f.seek(k)
#             dt = np.fromfile(f, dtype="float32", count=1)
#             k += num_bytes[4]
#             f.seek(k)
#             ns = np.fromfile(f, dtype="int32", count=1)[0]
#             k += num_bytes[5]
#             f.seek(k)
#             pcal = np.fromfile(f, dtype="float32", count=1)
#             k += num_bytes[6]
#             f.seek(k)
#             tcal = np.fromfile(f, dtype="float32", count=1)
#             k += num_bytes[7]
#             f.seek(k)
#             sm = np.fromfile(f, dtype="float32", count=1)
#             k += num_bytes[8]
#             f.seek(k)

#             # Only interested in the first component (the vertical component)
#             if j == 0:
#                 data_matrix[i, i_r, :] = np.fromfile(f, dtype="float32", count=ns)

#             k += num_bytes[9]
#             k += num_bytes[10]

#     nt = n_rec
#     ncomp = n_comp
#     f.close()

# with h5py.File(filenameH5Output, "w") as hf:

#     # Create group and attributes
#     g = hf.create_group("Texas synthetic data")
#     g.attrs["number of receivers"] = nt
#     g.attrs["number of components"] = ncomp
#     hf.create_dataset("stike, dip, rake, M0", data=sou_physics, dtype="int32")
#     hf.create_dataset("source location", data=sou_coordinates, dtype="int32")
#     hf.create_dataset("cluster IDs", data=one_hot_vectors)
#     hf.create_dataset("ML array", data=data_matrix, dtype="f")

""" Reading into numpy array """
import matplotlib.pyplot as plt

ncomp = 0
# Read single file to get number of receivers to initialize data_matrix
f = open(join(erz_folder, erzFiles[0]), "rb")
k = 4
f.seek(k)
nt = np.fromfile(f, dtype="int32", count=1)[0]  # number of receivers
f.close()

f = open(join(erz_folder, erzFiles[0]), "rb")
cmd_file = join(cmd_folder, erzFiles[0].split(".")[0] + ".cmd")
f_cmd = open(cmd_file, "r")

lines = f_cmd.read().splitlines()
f_cmd.close()

ns = int(lines[12].split()[0])

data = np.zeros((nt, ns))

# First part information about number of receivers and components per receiver
k = 4
f.seek(k)
n_rec = np.fromfile(f, dtype="int32", count=1)[0]  # number of receivers
k += 4
f.seek(k)
n_comp = np.fromfile(f, dtype="int32", count=1)[0]  # number of components per receiver
k += 8

# Not best prgramming. But does the job. These are the bytes at which to read beta info and data:
num_bytes = np.array([4, 4, 5, 3, 4, 4, 4, 4, 4, ns * 4, 4])

# Loop over all the receivers and their individual components
dists = []
azis = []
comps = []

st_syn = obspy.Stream()
for i_r in range(0, n_rec):
    for j in range(0, n_comp):

        k += num_bytes[0]
        f.seek(k)
        dist = np.fromfile(f, dtype="float32", count=1)
        k += num_bytes[1]
        f.seek(k)
        azi = np.fromfile(f, dtype="float32", count=1)
        k += num_bytes[2]
        f.seek(k)
        comp = np.fromfile(f, dtype="|S1", count=1).astype(str)[0]
        comps.append(comp)
        k += num_bytes[3]
        f.seek(k)
        dt = np.fromfile(f, dtype="float32", count=1)
        k += num_bytes[4]
        f.seek(k)
        ns = np.fromfile(f, dtype="int32", count=1)[0]
        k += num_bytes[5]
        f.seek(k)
        pcal = np.fromfile(f, dtype="float32", count=1)
        k += num_bytes[6]
        f.seek(k)
        tcal = np.fromfile(f, dtype="float32", count=1)
        k += num_bytes[7]
        f.seek(k)
        sm = np.fromfile(f, dtype="float32", count=1)
        k += num_bytes[8]
        f.seek(k)

        stats = {
            "delta": dt,
            "channel": "XB" + comp,
            "npts": ns,
            "starttime": or_time,
            "distance": dist * 1000,
        }
        tr_syn = obspy.Trace(np.fromfile(f, dtype="float32", count=ns), stats)
        st_syn += tr_syn

        # Only interested in the first component (the vertical component)
        if j == 2:

            data[i_r, :] = np.fromfile(f, dtype="float32", count=ns)

        k += num_bytes[9]
        k += num_bytes[10]
    dists.append(dist)
    azis.append(azi)

# st_select = st_syn.select(channel="XB" + comp)
# st_select.plot(type="section")

# fig, ax = plt.subplots(nrows=n_comp, ncols=1, figsize=(28, 12))
# for i, comp in enumerate(["Z", "R", "T"]):
#     st_select = st_syn.select(channel="XB" + comp)
#     for j, tr in enumerate(st_select):

#         ax[i].plot(tr.times(), tr.data, label=f"dist: {dists[j]}, azi: {azis[j]} ")
# plt.legend()
# plt.show()


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(28, 12))

colours = ["r", "b", "k"]

for i, comp in enumerate(["Z", "R", "T"]):
    st_select = st_syn.select(channel="XB" + comp)
    for j, tr in enumerate(st_select):

        ax.plot(
            tr.times(), tr.data + j * 1e-2, colours[i], label=f"dist: {dists[j]}, azi: {azis[j]} "
        )
# plt.legend()
plt.show()

