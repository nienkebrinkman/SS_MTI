import erzsol3Py.erzsol3Py as erz

""" 1.  Create csv-model file """
import numpy as np

# depth_model = np.arange(0, 3450, 50)
# dz = (depth_model[1] - depth_model[0]) / 1000  # convert to km
# Vp_model = np.linspace(5, 5, len(depth_model))
# Vs_model = Vp_model / 2

dz = np.array([1.0])
Vp_model = np.array([5.0])
Vs_model = np.array([2.6])
rho_model = np.ones_like(Vp_model) * 2.7
qa = 0.001
qb = 0.002

# """ 2.  Generating a .mod model """
folder = "/home/nienke/Documents/Research/SS_MTI/External_packages/erzsol3/Test_Files/"
mod_file_name = folder + "Test.mod"
file = open(mod_file_name, "w+")

model_name = "Test_model\n"
nLayers_thickness = "    {}        {}\n".format(len(Vp_model), 0)
L = [model_name, nLayers_thickness]
file.writelines(L)

for i, (vp, vs, rho) in enumerate(zip(Vp_model, Vs_model, rho_model)):
    L = [
        "3   {:.3f}    {:.3f}     {:.2f}     {:.3f}    {:.3f}     {:.3f}\n".format(
            vp, vs, rho, dz[i], qa, qb
        )
    ]
    file.writelines(L)

file.close()

""" Write .dst file """
epis = [22.0]
azis = [45.0]

dst_file_name = folder + "Test.dst"
file = open(dst_file_name, "w")
n_rec = len(epis)
L = ["   {}                               # of distances /distances\n".format(n_rec)]
file.writelines(L)

for i, (epi, azi) in enumerate(zip(epis, azis)):
    L = ["  {:.2f}      {:.2f}\n".format(epi, azi)]
    file.writelines(L)

file.close()

