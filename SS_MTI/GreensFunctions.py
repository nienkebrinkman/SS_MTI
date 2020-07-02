import obspy
import instaseis
from typing import Union as _Union
import numpy as np


def make_GF(
    self,
    or_time: obspy.UTCDateTime,
    lat_src: float,
    lon_src: float,
    depth: float,
    distance: float,
    rec: instaseis.Receiver,
    db: instaseis.open_db,
    dt: float,
    comp: str,
    tstar: _Union[float, str] = None,
    LQT: bool = False,
    inc: float = None,
    baz: float = None,
) -> obspy.Stream:
    """
    Create stream of different source components
    :param or_time: origin time
    :param lat_src: source latitude
    :param lon_src: source longitude
    :param depth: depth of event in km
    :param distance: the epicentral distance in degrees
    :param rec: instaseis.Receiver object of the single station
    :param db: instaseis database
    :param dt: timestep
    :param comp: component
    :param tstar: tstar value 
    :param LQT: set to true if component system is LQT
    :param inc: inclination angle in degrees (needed when LQT = TRUE)
    :param baz: backazimuth angle in degrees (needed when LQT = TRUE)
    """

    if tstar is not None and not isinstance(tstar, str):
        stf_len_sec = 30.0
        stf = self.STF_.stf_tstar(tstar=tstar, dt=db.info.dt, npts=int(stf_len_sec / db.info.dt))[
            0
        ]
    elif isinstance(tstar, str):
        stf = self.STF_.Create_stf_from_file(tstar, db.info.dt)
    mts = [
        [1e14, 0e14, 0e14, 0e14, 0e14, 0e14],
        [0e14, 1e14, 0e14, 0e14, 0e14, 0e14],
        [0e14, 0e14, 1e14, 0e14, 0e14, 0e14],
        [0e14, 0e14, 0e14, 1e14, 0e14, 0e14],
        [0e14, 0e14, 0e14, 0e14, 1e14, 0e14],
        [0e14, 0e14, 0e14, 0e14, 0e14, 1e14],
    ]

    st = obspy.Stream()

    for mt in mts:
        src = instaseis.Source(
            latitude=lat_src,
            longitude=lon_src,
            depth_in_m=depth * 1e3,
            origin_time=or_time,
            m_rr=mt[0],
            m_pp=mt[1],
            m_tt=mt[2],
            m_rp=mt[3],
            m_rt=mt[4],
            m_tp=mt[5],
        )

        reconvolve_stf = False
        remove_source_shift = True
        if tstar is not None and not isinstance(tstar, str):
            reconvolve_stf = True
            remove_source_shift = False
            src.set_sliprate(stf, dt=db.info.dt)
            reconvolve_stf = True
            remove_source_shift = False
        elif isinstance(tstar, str):
            reconvolve_stf = True
            remove_source_shift = False
            src.set_sliprate(stf, dt=db.info.dt, normalize=True)

        if LQT:
            st_rot = db.get_seismograms(
                src,
                rec,
                dt=dt,
                components="ZNE",
                kind="displacement",
                reconvolve_stf=reconvolve_stf,
                remove_source_shift=remove_source_shift,
            )
            st_rot.rotate(method="ZNE->LQT", back_azimuth=baz, inclination=inc)
            tr_rot = st_rot.select(channel="BX" + comp[0])[0]
            st += tr_rot
        else:
            st += db.get_seismograms(
                src,
                rec,
                dt=dt,
                components=comp,
                kind="displacement",
                reconvolve_stf=reconvolve_stf,
                remove_source_shift=remove_source_shift,
            )[0]
    return st


def convert_SDR(strike, dip, rake, M0=1e14):
    phi = np.deg2rad(strike)
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    m_rr = (np.sin(2.0 * delta) * np.sin(lambd)) * M0

    m_pp = (
        np.sin(delta) * np.cos(lambd) * np.sin(2.0 * phi)
        - np.sin(2.0 * delta) * np.cos(phi) ** 2.0 * np.sin(lambd)
    ) * M0

    m_tt = (
        -np.sin(delta) * np.cos(lambd) * np.sin(2.0 * phi)
        - np.sin(2.0 * delta) * np.sin(phi) ** 2.0 * np.sin(lambd)
    ) * M0

    m_rp = (
        -np.cos(phi) * np.sin(lambd) * np.cos(2.0 * delta)
        + np.cos(delta) * np.cos(lambd) * np.sin(phi)
    ) * M0

    m_rt = (
        -np.sin(lambd) * np.sin(phi) * np.cos(2.0 * delta)
        - np.cos(delta) * np.cos(lambd) * np.cos(phi)
    ) * M0

    m_tp = (
        -np.sin(delta) * np.cos(lambd) * np.cos(2.0 * phi)
        - np.sin(2.0 * delta) * np.sin(2.0 * phi) * np.sin(lambd) / 2.0
    ) * M0

    MT = [m_rr, m_pp, m_tt, m_rp, m_rt, m_tp]
    return MT


def from_GF(st_in, strike, dip, rake, M0=1e14):
    MT = convert_SDR(strike, dip, rake, M0)
    m_rr = MT[0] / 1e14
    m_pp = MT[1] / 1e14
    m_tt = MT[2] / 1e14
    m_rp = MT[3] / 1e14
    m_rt = MT[4] / 1e14
    m_tp = MT[5] / 1e14

    data = (
        st_in[0].data * m_rr
        + st_in[1].data * m_pp
        + st_in[2].data * m_tt
        + st_in[3].data * m_rp
        + st_in[4].data * m_rt
        + st_in[5].data * m_tp
    )

    tr = st_in[0].copy()
    tr.data = data

    return tr
