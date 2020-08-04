import SS_MTI
import obspy
import EventInterface
import instaseis

## Parameters:
or_time = obspy.UTCDateTime("2020-3-10T12:00:00")
lat_src = 10.99032013
lon_src = 170
depth = 45.0
name = "Test_Event"

phases = ["P", "S"]
lat_rec = 4.502384
lon_rec = 135.623447
rec = instaseis.Receiver(
    latitude=lat_rec, longitude=lon_rec, network="XB", station="ELYSE", location="02",
)
epi, az, baz = EventInterface.EventObj.Get_location(
    la_s=lat_src, lo_s=lon_src, la_r=lat_rec, lo_r=lon_rec
)

mnt_folder = "/mnt/marshost/"

SS_MTI.DataGetter.mnt_remote_folder(
    host_ip="marshost.ethz.ch", host_usr="sysop", remote_folder="/data/", mnt_folder=mnt_folder,
)

db_path = "/mnt/marshost/instaseis2/databases/TAYAK_15s_BKE"
# db_path = "http://instaseis.ethz.ch/blindtest_1s/TAYAK_1s/"
db = instaseis.open_db(db_path)

dt = 0.05


def Get_GF_with_STF(origin_time, baz, tstar, db, epi, depth_in_m, dt, LQT_value=False, inc=None):
    src_latitude, src_longitude = 90.0, 0.0
    rec_latitude, rec_longitude = 90.0 - epi, 0.0

    # sources according to https://github.com/krischer/instaseis/issues/8
    # transformed to r, theta, phi
    #
    # Mtt =  Mxx, Mpp = Myy, Mrr =  Mzz
    # Mrp = -Myz, Mrt = Mxz, Mtp = -Mxy
    #
    # Mrr   Mtt   Mpp    Mrt    Mrp    Mtp
    #  0     0     0      0      0     -1.0    m1
    #  0     1.0  -1.0    0      0      0      m2
    #  0     0     0      0     -1.0    0      m3
    #  0     0     0      1.0    0      0      m4
    #  1.0   1.0   1.0    0      0      0      m6
    #  2.0  -1.0  -1.0    0      0      0      cl

    m1 = instaseis.Source(
        src_latitude, src_longitude, depth_in_m, m_tp=-1.0, origin_time=origin_time
    )
    m2 = instaseis.Source(
        src_latitude, src_longitude, depth_in_m, m_tt=1.0, m_pp=-1.0, origin_time=origin_time
    )
    m3 = instaseis.Source(
        src_latitude, src_longitude, depth_in_m, m_rp=-1.0, origin_time=origin_time
    )
    m4 = instaseis.Source(
        src_latitude, src_longitude, depth_in_m, m_rt=1.0, origin_time=origin_time
    )
    m6 = instaseis.Source(
        src_latitude,
        src_longitude,
        depth_in_m,
        m_rr=1.0,
        m_tt=1.0,
        m_pp=1.0,
        origin_time=origin_time,
    )
    cl = instaseis.Source(
        src_latitude,
        src_longitude,
        depth_in_m,
        m_rr=2.0,
        m_tt=-1.0,
        m_pp=-1.0,
        origin_time=origin_time,
    )

    receiver = instaseis.Receiver(rec_latitude, rec_longitude)

    if LQT_value:
        items = [
            ("SST", m1, "T"),
            ("SSL", m2, "L"),
            ("SSQ", m2, "Q"),
            ("DST", m3, "T"),
            ("DSL", m4, "L"),
            ("DSQ", m4, "Q"),
            ("DDL", cl, "L"),
            ("DDQ", cl, "Q"),
            ("EPL", m6, "L"),
            ("EPQ", m6, "Q"),
        ]

    else:
        items = [
            ("SST", m1, "T"),
            ("SSZ", m2, "Z"),
            ("SSR", m2, "R"),
            ("DST", m3, "T"),
            ("DSZ", m4, "Z"),
            ("DSR", m4, "R"),
            ("DDZ", cl, "Z"),
            ("DDR", cl, "R"),
            ("EPZ", m6, "Z"),
            ("EPR", m6, "R"),
        ]

    args = {"receiver": receiver, "dt": dt, "kind": "displacement", "kernelwidth": 12}
    st = obspy.Stream()

    if tstar is not None and not isinstance(tstar, str):
        stf_len_sec = 30.0
        stf = self.STF_.stf_tstar(tstar=tstar, dt=db.info.dt, npts=int(stf_len_sec / db.info.dt))[
            0
        ]
    elif isinstance(tstar, str):
        stf = self.STF_.Create_stf_from_file(tstar, db.info.dt)

    for name, src, comp in items:
        reconvolve_stf = False
        remove_source_shift = True
        if tstar is not None and not isinstance(tstar, str):
            src.set_sliprate(stf, dt=db.info.dt)
            reconvolve_stf = True
            remove_source_shift = False
        elif isinstance(tstar, str):
            src.set_sliprate(stf, dt=db.info.dt, normalize=True)
            reconvolve_stf = True
            remove_source_shift = False

        if LQT_value:
            st_rot = db.get_seismograms(
                source=src,
                components="ZRT",
                reconvolve_stf=reconvolve_stf,
                remove_source_shift=remove_source_shift,
                **args
            )
            st_rot.rotate(method="RT->NE", back_azimuth=baz)
            st_rot.rotate(method="ZNE->LQT", back_azimuth=baz, inclination=inc)
            tr = st_rot.select(channel="BX" + comp)[0]
        else:
            tr = db.get_seismograms(
                source=src,
                components=comp,
                reconvolve_stf=reconvolve_stf,
                remove_source_shift=remove_source_shift,
                **args
            )[0]

        tr.stats.channel = name
        st.append(tr)
    return st


st = SS_MTI.GreensFunctions.make_GF(
    or_time=or_time,
    lat_src=lat_src,
    lon_src=lon_src,
    depth=depth,
    distance=epi,
    rec=rec,
    db=db,
    dt=dt,
    comp="T",
    tstar=None,
    LQT=False,
    inc=None,
    baz=None,
    M0=1e14,
)

st_test = Get_GF_with_STF(
    origin_time=or_time,
    baz=baz,
    tstar=None,
    db=db,
    epi=epi,
    depth_in_m=depth * 1e3,
    dt=dt,
    LQT_value=False,
    inc=None,
)

rec_in = instaseis.Receiver(
    latitude=90.0 - epi, longitude=0.0, network="XB", station="ELYSE", location="02",
)
st_in = SS_MTI.GreensFunctions.make_GF(
    or_time=or_time,
    lat_src=90.0,
    lon_src=0.0,
    depth=depth,
    distance=epi,
    rec=rec_in,
    db=db,
    dt=dt,
    comp="T",
    tstar=None,
    LQT=False,
    inc=None,
    baz=None,
    M0=1.0,
)
SS, DS = SS_MTI.GreensFunctions.from_GF_get_G(st_in=st_in, az=az, comp="T")
import matplotlib.pyplot as plt

plt.plot(st_test.traces[0].data)
plt.plot(SS)
plt.show()

a = 1

