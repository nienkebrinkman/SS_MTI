from obspy.taup import TauPyModel as TauPyModel
from os.path import isfile, join
import numpy as np


def create_dat_file(
    src_depth: float,
    epi: float,
    baz: float,
    focal_mech: [float],
    M0: float = None,
    save_path: str = "./",
    bm_file_path: str = "./",
    fdom: str = 1.000,
):
    """ 
    This function creates a .dat file that is used for the reflectivity code of Fuchs&Muller
    :paran src_depth: source depth
    :param focal_mech: strike,dip,rake or m_rr, m_tt, m_pp, m_rt, m_rp, m_tp
    :param M0: scalar moment, only necessesary when focal_mech strike,dip,rake
    :param epi: Epicentral distance (in degrees)
    :param baz: Back-azimuth (in degrees)
    :param save_path: path to save .dat file
    :param fdom: dominant frequency
    """

    bm_file = bm_file_path

    f = np.loadtxt(bm_file, skiprows=5)
    f_ud = np.flipud(f)

    radius_mars = f_ud[0][0]  # 3390 (km)

    # radius_of_planet = 3390
    km_per_deg = np.pi * (radius_mars * 1e-3) / 180.0
    dist_in_km = epi * np.pi * (radius_mars * 1e-3) / 180.0
    dist = dist_in_km

    if baz < 0:
        baz *= -1
    rec_az = baz
    rec_z = 0.0

    src_x = 0.0
    src_y = 0.0
    src_z = src_depth
    or_time = 0.0
    s_strength = 1.0

    assert (M0 is None and len(focal_mech) == 6) or (M0 is not None and len(focal_mech) == 3), (
        "focal_mech length is incorrect. "
        "If you specify M0, focal_mech is [strike,dip,rake]. "
        "Otherwise focal_mech is [m_rr, m_tt, m_pp, m_rt, m_rp, m_tp]"
    )

    for i in range(len(focal_mech)):
        focal_mech[i] += 0

    M_tt_ins = focal_mech[1]
    M_pp_ins = focal_mech[2]
    M_rr_ins = focal_mech[0]
    M_rp_ins = focal_mech[4]
    M_rt_ins = focal_mech[3]
    M_tp_ins = focal_mech[5]

    moment_tensor = f"{M_tt_ins:10.4f}{M_tp_ins:10.4f}{-M_rt_ins+0:10.4f}{M_pp_ins:10.4f}{-M_rp_ins+0:10.4f}{M_rr_ins:10.4f}"
    # moment_tensor = f"{M_tt_ins:10.4f}{M_tp_ins:10.4f}{M_rt_ins:10.4f}{M_pp_ins:10.4f}{M_rp_ins:10.4f}{M_rr_ins:10.4f}"

    # model = TauPyModel(taup_path)
    # model_layers = model.model.s_mod.v_mod.layers

    with open(join(save_path, "crfl.dat"), "w") as f:
        f.write("Test name\n")
        f.write(" 0 0 0 0 0   0 0 1 1 1   2 1 0 0 1   0 1 2 0 1   1\n")
        f.write("    5    1    0    1    1\n")

        for layer in f_ud:
            if layer[0] == 0.0:
                continue
            depth = (radius_mars - layer[0]) * 1e-3
            dens = layer[1] * 1e-3
            vp = layer[2] * 1e-3
            vs = layer[3] * 1e-3
            qka = layer[4]  # qka
            qmu = layer[5]  # qmu
            vph = layer[6]
            vsh = layer[7]
            eta = layer[8]

            qs = qmu
            L = (4 / 3) * (vs / vp) ** 2
            qp = 1 / (L * (1 / qmu) + (1 - L) * (1 / qka))
            if np.isnan(qp):
                qp = qka
                qs = 10.0
            text = f"{depth:10.4f}{vp:10.4f}{qp:10.4f}{vs:10.4f}{qs:10.4f}{dens:10.4f}{1:10d}\n"
            f.write(text)

        #     with open(join(save_path, "crfl.dat"), "w") as f:
        #         f.write("Test name\n")
        #         f.write(" 0 0 0 0 0   0 0 1 1 1   2 1 0 0 1   0 1 2 0 1   1\n")
        #         f.write("    5    1    0    1    1\n")
        #         f.write(
        #             """    0.0000    3.6777 1499.3884    1.7398  500.0000    1.8653        30
        #     1.0000    3.6777 1499.4177    1.7398  500.0000    1.8653         1
        #     1.0000    4.9523 1112.4995    2.7810  500.0000    2.2718        30
        #    10.0000    4.9523 1112.7626    2.7806  500.0000    2.2718         1
        #    10.0000    5.8467 1113.7485    3.2812  500.0000    2.6817        30
        #    77.3680    5.8467 1115.6692    3.2780  500.0000    2.6817         1
        #    77.3680    7.4009  265.2819    4.2463  118.2000    3.3886         1
        #    80.0000    7.4009  265.2819    4.2463  118.2000    3.3886         1
        #    80.0000    7.4136  266.0473    4.2473  118.2000    3.3931         1
        #   100.0000    7.4361  276.5086    4.1802  118.4000    3.4000         1
        #   110.0000    7.4466  277.1442    4.1812  118.4000    3.4033         1
        #   120.0000    7.4985  279.2070    4.1960  118.5000    3.4154         1
        #   130.0000    7.5241  280.6065    4.2013  118.6000    3.4216         1
        #   140.0000    7.5466  282.1372    4.2058  118.8000    3.4273         1
        #   150.0000    7.5647  283.2415    4.2091  118.9000    3.4321         1
        #   160.0000    7.5833  284.3444    4.2128  119.0000    3.4366         1
        #   170.0000    7.6026  285.6435    4.2173  119.2000    3.4414         1
        #   180.0000    7.6158  286.5863    4.2192  119.3000    3.4452         1
        #   190.0000    7.6281  287.6835    4.2213  119.5000    3.4489         1
        #   212.6800    7.6609  290.3808    4.2263  119.9000    3.4578         1
        #   235.0760    7.6856  292.6324    4.2302  120.3000    3.4654         1
        #   257.1030    7.7110  295.2375    4.2337  120.8000    3.4732         1
        #   279.1150    7.7323  297.4730    4.2360  121.2000    3.4800         1
        #   301.1010    7.7542  300.1527    4.2390  121.8000    3.4869         1
        #   323.1030    7.7782  302.6655    4.2427  122.3000    3.4944         1
        #   345.1490    7.7974  305.2588    4.2450  122.9000    3.5007         1
        #   366.8120    7.8182  307.9317    4.2476  123.5000    3.5073         1
        #   386.4420    7.8338  310.3923    4.2491  124.1000    3.5129         1
        #   406.4120    7.8507  312.8767    4.2512  124.7000    3.5185         1
        #   426.3030    7.8682  315.4224    4.2532  125.3000    3.5245         1
        #   445.8490    7.8836  317.8886    4.2547  125.9000    3.5299         1
        #   465.4600    7.8992  320.3666    4.2563  126.5000    3.5354         1
        #   485.0010    7.9144  323.1287    4.2574  127.2000    3.5409         1
        #   504.5090    7.9368  319.8783    4.2662  125.7000    3.5477         1
        #   524.0180    7.9575  317.9313    4.2735  124.7000    3.5542         1
        #   543.5260    7.9778  316.2334    4.2805  123.8000    3.5606         1
        #   563.0340    8.0008  314.8885    4.2882  123.0000    3.5678         1
        #   582.5430    8.0208  313.3795    4.2954  122.2000    3.5741         1
        #   602.0510    8.0408  312.1230    4.3025  121.5000    3.5804         1
        #   621.5600    8.0610  310.8379    4.3099  120.8000    3.5868         1
        #   641.0680    8.0802  309.8039    4.3167  120.2000    3.5931         1
        #   660.5770    8.0988  309.0361    4.3231  119.7000    3.5992         1
        #   680.0850    8.1185  308.5141    4.3300  119.3000    3.6055         1
        #   699.5930    8.1378  307.9662    4.3369  118.9000    3.6118         1
        #   719.1020    8.1579  307.6540    4.3444  118.6000    3.6181         1
        #   738.6100    8.1769  307.5782    4.3513  118.4000    3.6243         1
        #   758.1190    8.1965  307.5002    4.3586  118.2000    3.6308         1
        #   777.6270    8.2155  307.6610    4.3656  118.1000    3.6369         1
        #   797.1360    8.2346  307.8415    4.3726  118.0000    3.6432         1
        #   816.6440    8.2548  308.2487    4.3803  118.0000    3.6498         1
        #   836.1530    8.2738  308.9393    4.3872  118.1000    3.6562         1
        #   855.6610    8.2940  309.5893    4.3950  118.2000    3.6628         1
        #   875.1690    8.3453  308.7424    4.4304  118.3000    3.6724         1
        #   894.6780    8.3649  309.6833    4.4376  118.5000    3.6794         1
        #   914.1860    8.3834  310.6247    4.4442  118.7000    3.6859         1
        #   933.6950    8.4041  311.5517    4.4522  118.9000    3.6935         1
        #   953.2030    8.4228  313.0306    4.4588  119.3000    3.7005         1
        #   972.7120    8.4443  314.5816    4.4663  119.7000    3.7087         1
        #   992.2200    8.4658  316.1571    4.4737  120.1000    3.7170         1
        #  1011.7280    8.4862  318.1342    4.4795  120.6000    3.7257         1
        #  1031.2370    8.5863  317.3462    4.5381  120.6000    3.7550         1
        #  1050.7450    8.7199  315.4967    4.6168  120.3000    3.7935         1
        #  1070.2540    8.9089  312.3254    4.7336  119.9000    3.8441         1
        #  1089.7620    9.0000  311.6404    4.7874  119.9000    3.8686         1
        #  1109.2710    9.0225  313.2845    4.7965  120.4000    3.8767         1
        #  1128.7790    9.0400  314.8755    4.8032  120.9000    3.8828         1
        #  1148.2880    9.0585  316.7598    4.8103  121.5000    3.8899         1
        #  1167.7960    9.0745  318.6537    4.8160  122.1000    3.8957         1
        #  1187.3040    9.0896  320.7260    4.8218  122.8000    3.9010         1
        #  1206.8130    9.1073  323.0651    4.8290  123.6000    3.9075         1
        #  1226.3210    9.1225  325.3992    4.8348  124.4000    3.9131         1
        #  1245.8300    9.1607  327.0772    4.8540  125.0000    3.9227         1
        #  1265.3380    9.2065  328.1739    4.8777  125.4000    3.9346         1
        #  1284.8470    9.2408  329.5776    4.8950  125.9000    3.9436         1
        #  1304.3550    9.2873  330.5675    4.9199  126.3000    3.9558         1
        #  1323.8630    9.3329  331.5513    4.9443  126.7000    3.9676         1
        #  1343.3720    9.3830  332.5809    4.9728  127.2000    3.9808         1
        #  1362.8800    9.4277  333.0273    4.9970  127.4000    3.9935         1
        #  1382.3890    9.4643  334.9041    5.0177  128.2000    4.0039         1
        #  1401.8970    9.4796  338.3990    5.0245  129.5000    4.0095         1
        #  1421.4060    9.4946  341.8790    5.0313  130.8000    4.0149         1
        #  1440.9140    9.5091  345.6042    5.0379  132.2000    4.0201         1
        #  1460.4230    9.5237  349.8997    5.0441  133.8000    4.0251         1
        #  1479.9310    9.5381  354.1162    5.0508  135.4000    4.0301         1
        #  1499.4390    9.5511  358.6519    5.0563  137.1000    4.0351         1
        #  1518.9480    9.5658  363.4052    5.0630  138.9000    4.0408         1
        #  1538.4560    9.5801  368.6329    5.0697  140.9000    4.0458         1
        #  1557.9650    9.5943  374.1197    5.0764  143.0000    4.0508         1
        #  1577.4730    9.6084  380.1276    5.0829  145.3000    4.0558         1
        #  1596.9820    9.6379  384.4757    5.1020  147.2000    4.0668         1
        #  1596.9820    4.986810000.0000    0.0000   10.0000    5.8554         1
        #  1853.0560    5.150010000.0000    0.0000   10.0000    6.0257         1
        #  2109.1300    5.284810000.0000    0.0000   10.0000    6.1669         1
        #  2365.2040    5.393010000.0000    0.0000   10.0000    6.2807         1
        #  2621.2780    5.475910000.0000    0.0000   10.0000    6.3682         1
        #  2877.3520    5.534510000.0000    0.0000   10.0000    6.4301         1
        #  3133.4260    5.569510000.0000    0.0000   10.0000    6.4671         1
        # """
        #         )

        f.write("\n")
        f.write(f"{rec_z:10.4f}\n")
        f.write(f"{src_x:10.4f}{src_y:10.4f}{src_z:10.4f}{or_time:10.4f}{s_strength:10.4f}\n")
        f.write(f"{moment_tensor}\n")
        f.write(f"{dist:10.4f}{dist:10.4f}{0.:10.4f}{rec_az:10.4f}{1:10d}\n")
        f.write(f"{dist:10.4f}\n")
        f.write(f"{rec_az:10.4f}\n")
        f.write(f"{12.:10.4f}   {-300.:10.4f}\n")
        f.write("    3.0000    3.5000   23.5000   25.0000      650\n")
        f.write("    0.0100    0.0133    1.0000    1.0300    0.0000\n")
        f.write("    0.0250     65536         0         2    0.0250  491.5200\n")

    f.close()


save_path = "/home/nienke/Documents/Research/SS_MTI/External_packages/Test_reflectivity/"

bm_file_path = "/home/nienke/Documents/Research/Data/MTI/MT_vs_STR/bm_models/TAYAK.bm"
create_dat_file(
    src_depth=40.0,
    focal_mech=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    M0=None,
    epi=30.0,
    baz=0.0,
    save_path=save_path,
    bm_file_path=bm_file_path,
)
