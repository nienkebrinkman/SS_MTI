#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    Nienke Brinkman (nienke.brinkman@erdw.ethz.ch), 2020
:license:
    None
"""

import obspy
import instaseis
from typing import Union as _Union
import numpy as np

import SS_MTI.SourceTimeFunction as _STF


def make_GF(
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
    M0: float = 1e14,
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
    :param M0: scalar moment
    """

    if tstar is not None and not isinstance(tstar, str):
        stf_len_sec = 30.0
        stf = _STF.stf_tstar(tstar=tstar, dt=db.info.dt, npts=int(stf_len_sec / db.info.dt))[0]
    elif isinstance(tstar, str):
        stf = _STF.Create_stf_from_file(tstar, db.info.dt)
    mts = [
        [M0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, M0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, M0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, M0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, M0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, M0],
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


def convert_SDR(strike: float, dip: float, rake: float, M0: float = 1e14):
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


def from_GF(st_in: obspy.Stream, focal_mech: [float], M0: float):
    """ Generate synthetic waveforms 
    :param st_in: 
    :param focal_mech: strike,dip,rake or m_rr, m_pp, m_tt, m_rp, m_rt, m_tp
    :param M0: scalar moment
    """

    if len(focal_mech) == 3:
        focal_mech = convert_SDR(focal_mech[0], focal_mech[1], focal_mech[2], M0)

    m_rr = focal_mech[0] / M0
    m_pp = focal_mech[1] / M0
    m_tt = focal_mech[2] / M0
    m_rp = focal_mech[3] / M0
    m_rt = focal_mech[4] / M0
    m_tp = focal_mech[5] / M0

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


def from_GF_get_G(st_in: obspy.Stream, az: float, comp: str):
    m_rr = st_in[0].data
    m_pp = st_in[1].data
    m_tt = st_in[2].data
    m_rp = st_in[3].data
    m_rt = st_in[4].data
    m_tp = st_in[5].data

    m1 = -1.0 * m_tp
    m2 = 1.0 * m_tt + -1.0 * m_pp
    m3 = -1.0 * m_rp
    m4 = 1.0 * m_rt
    m6 = 1.0 * m_rr + 1.0 * m_tt + 1.0 * m_pp
    cl = 2.0 * m_rr + -1.0 * m_tt + -1.0 * m_pp

    if comp == "Z" or comp == "R" or comp == "L" or comp == "Q":
        SS = m2
        DS = m4
        DD = cl
        EP = m6  # Explosion term

        G = np.zeros((len(SS), 5))
        G[:, 0] = SS * (0.5) * np.cos(2 * np.deg2rad(az)) - DD / 2.0
        G[:, 1] = -SS * (0.5) * np.cos(2 * np.deg2rad(az)) - DD / 2.0
        G[:, 2] = SS * np.sin(2 * np.deg2rad(az))
        G[:, 3] = -DS * np.cos(np.deg2rad(az))
        G[:, 4] = -DS * np.sin(np.deg2rad(az))

    elif comp == "T":
        SS = m1
        DS = m3

        G = np.zeros((len(SS), 5))
        G[:, 0] = -SS * (0.5) * np.sin(2 * np.deg2rad(az))
        G[:, 1] = SS * (0.5) * np.sin(2 * np.deg2rad(az))
        G[:, 2] = SS * np.cos(2 * np.deg2rad(az))
        G[:, 3] = DS * np.sin(np.deg2rad(az))
        G[:, 4] = -DS * np.cos(np.deg2rad(az))
    else:
        raise ValueError("Component is not correctly specified")
    return G


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
        stf = _STF.stf_tstar(tstar=tstar, dt=db.info.dt, npts=int(stf_len_sec / db.info.dt))[0]
    elif isinstance(tstar, str):
        stf = _STF.Create_stf_from_file(tstar, db.info.dt)

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

