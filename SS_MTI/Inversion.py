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
import numpy as _np
import h5py as _h5
import matplotlib.pyplot as plt
from os.path import join as pjoin
from obspy.taup import TauPyModel
from obspy import UTCDateTime as utct
from typing import List as _List, Union as _Union


from SS_MTI import PreProcess as _PreProcess
from SS_MTI import Forward as _Forward
from SS_MTI import Misfit as _Misfit
from SS_MTI import MTDecompose as _MTDecompose
from SS_MTI import PostProcessing as _PostProcessing


def Grid_Search_run(
    fwd: _Forward._AbstractForward,
    misfit: _Misfit._AbstractMisfit,
    event: obspy.core.event.Event,
    rec: instaseis.Receiver,
    phases: [str],
    components: [str],
    t_pre: [str],
    t_post: [str],
    depths: [float],
    strikes: [float],
    dips: [float],
    rakes: [float],
    phase_corrs: [float] = None,
    tstars: _Union[_List[float], _List[str]] = None,
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = False,
    list_to_correct_M0: [str] = None,
    output_folder: str = None,
    plot: bool = False,
    plot_extra_phases: [str] = None,
    color_plot: str = None,
    Ylims: [float] = None,
):
    """
    Grid search over strike, dip, rake angles
    :param event: Obspy.event including waveforms and phase arrivals
    :param phases: list of phases to include in the inversion
    """
    print(f"Running grid search with model: {fwd.name}")
    print(f"and with {misfit.description}")
    M0 = 1.0

    if tstars is None:
        tstars = [None] * len(phases)

    if phase_corrs is None:
        phase_corrs = [0] * len(phases)

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    if output_folder is None:
        output_folder = "."

    # TODO: IMPLEMENT LQT COORDINATE SYSTEM
    LQT_value = False
    baz = None
    inc = None

    """ PRE-PROCESS THE OBSERVED DATA """
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
    st_obs, sigmas = _PreProcess.prepare_event_data(
        event=event,
        phases=phases,
        components=components,
        slice=True,
        tts=obs_tt,
        t_pre=t_pre,
        t_post=t_post,
        filter=filter_par,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        noise_level=misfit.noise_level,
    )

    for depth in depths:
        print(depth)
        """ Open .h5 file """
        file_name = f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.hdf5"
        f = _h5.File(pjoin(output_folder, file_name), "w")
        data_len = 6 + len(phases)
        file_len = len(strikes) * len(dips) * len(rakes)
        f.create_dataset("samples", (file_len, data_len), maxshape=(None, 50))
        iteration = 0

        """ GENERATE GREEN'S FUNCTION AT SPECIFIC DEPTH """
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
                M0=M0,
                filter=filter_par,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
            )
            syn_GFs.append(syn_GF)
            syn_tts.append(syn_tt)

        for strike in strikes:
            for dip in dips:
                for rake in rakes:

                    focal_mech = [strike, dip, rake]
                    st_syn = obspy.Stream()
                    misfit_amp = []

                    """ Generate the synthetic data"""
                    for i, phase in enumerate(phases):
                        tr_syn = fwd.generate_synthetic_data(
                            st_GF=syn_GFs[i],
                            focal_mech=focal_mech,
                            M0=M0,
                            slice=True,
                            tt=syn_tts[i],
                            t_pre=t_pre[i],
                            t_post=t_post[i],
                        )

                        if phases[i] + components[i] in list_to_correct_M0:
                            misfit_amp.append(
                                (_np.sum(_np.abs(st_obs[i].data)))
                                / (_np.sum(_np.abs(tr_syn.data)))
                            )

                        st_syn += tr_syn

                    """ Multiply the data with the M0 correction"""
                    M0_corr = _np.sum(misfit_amp) / len(misfit_amp)
                    for tr in st_syn:
                        tr.data = tr.data * M0_corr

                    """ Determine the misfit between syntetic and observed"""
                    chi = misfit.run_misfit(
                        phases=phases, st_obs=st_obs, st_syn=st_syn, sigmas=sigmas
                    )

                    """ Write into file"""
                    f["samples"][iteration, :] = [depth, strike, dip, rake, M0, M0_corr] + chi
                    iteration += 1

                    # print(focal_mech)
                    # print(_np.sum(chi))
                    # print(M0_corr * M0)

        if plot:
            # TODO: make this plot routine as a function
            """ Calculate take-off angles"""
            takeoff_angles = ["P", "S", "pP"]
            angles = []
            for phase in takeoff_angles:
                angles.append(
                    fwd.get_phase_tt(
                        phase=phase, depth=depth, distance=event.distance, takeoffs=True
                    )
                )
            if color_plot is None:
                color_plot = "blue"

            sum_misfits = _np.sum(f["samples"][:, -len(phases) :], axis=1)
            nlowest = 1
            lowest_indices = sum_misfits.argsort()[0:nlowest]
            sdrs_total = f["samples"][:, 1:4]
            sdrs = sdrs_total[lowest_indices, :]
            M0_corrs_total = f["samples"][:, 5]
            M0_corrs = M0_corrs_total[lowest_indices]
            M0_total = f["samples"][:, 4]
            M0_plot = _np.expand_dims(
                M0_total[lowest_indices] * M0_corrs_total[lowest_indices], axis=1
            )
            depth_post = f["samples"][:, 0]
            depth_plot = depth_post[lowest_indices]
            """ Beachball plot """

            fig = _PostProcessing.Plot_GS_BB(
                sdrs[:, 0],
                sdrs[:, 1],
                sdrs[:, 2],
                azimuths=[event.az, event.az, event.az],
                inc_angles=angles,
                phase_names=takeoff_angles,
                color=color_plot,
            )
            plt.savefig(
                pjoin(
                    output_folder,
                    f"GS_BBB__{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}.svg",
                ),
                dpi=300,
            )
            plt.close()

            extra_arrs = []
            for j, extraphase in enumerate(plot_extra_phases):
                arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
                extra_arrs.append(arr)

            fig, ax = _PostProcessing.waveform_plot(
                syn_GFs=syn_GFs,
                syn_tts=syn_tts,
                obs_tts=obs_tt,
                fwd=fwd,
                misfit_weight_len=misfit.start_weight_len,
                event=event,
                phases=phases,
                components=components,
                t_pre=t_pre,
                t_post=t_post,
                MTs=sdrs,
                M0s=M0_plot,
                fmin=fmin,
                fmax=fmax,
                zerophase=zerophase,
                plot_extra_phases=plot_extra_phases,
                extra_arrs=extra_arrs,
                color_plot=color_plot,
                Ylims=Ylims,
            )

            for i, phase in enumerate(phases):
                for n in range(len(M0_plot)):
                    st = make_GF(
                        event=event,
                        depth=depth_plot[n],
                        rec=rec,
                        db=fwd.db,
                        dt=fwd.dt,
                        tstar=tstars[i],
                        comp=[components[i]],
                        LQT=False,
                        inc=None,
                    )
                    st.trim(
                        starttime=fwd.or_time + fwd.start_cut, endtime=fwd.or_time + fwd.end_cut
                    )
                    _PreProcess.filter_tr(st, fmin=fmin, fmax=fmax, zerophase=zerophase)

                    tr_syn_full = from_GF(
                        st_in=st, strike=sdrs[n, 0], dip=sdrs[n, 1], rake=sdrs[n, 2], M0=1.0,
                    )

                    tr_slice = tr_syn_full.slice(
                        starttime=fwd.or_time + syn_tts[i] - t_pre[i],
                        endtime=fwd.or_time + syn_tts[i] + t_post[i],
                    )

                    ax[i].plot(
                        tr_slice.times() - t_pre[i], tr_slice.data * M0_corrs[n], lw=0.5, c="r",
                    )
                    ax[i].plot(
                        tr_syn_full.times() - (syn_tts[i] - fwd.start_cut),
                        tr_syn_full.data * M0_corrs[n],
                        lw=0.5,
                        c="r",
                    )
            plt.savefig(
                pjoin(
                    output_folder,
                    f"GS_waveforms_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}.svg",
                ),
                dpi=300,
            )
            plt.close()
        f.close()


def make_GF(event, depth, rec, db, dt, tstar, comp, LQT=False, inc=None):

    """  Put tstar = None if you do not want to use a STF"""
    import SS_MTI.SourceTimeFunction as _STF

    if tstar is not None and not isinstance(tstar, str):
        stf_len_sec = 30.0
        stf = _STF.stf_tstar(tstar=tstar, dt=db.info.dt, npts=int(stf_len_sec / db.info.dt))[0]
    elif isinstance(tstar, str):
        stf = _STF.Create_stf_from_file(tstar, db.info.dt)

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
            latitude=event.latitude,
            longitude=event.longitude,
            depth_in_m=depth * 1e3,
            origin_time=event.origin_time,
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
            st_rot.rotate(method="ZNE->LQT", back_azimuth=event.baz, inclination=inc)
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


def from_GF(st_in, strike, dip, rake, M0=1e14):
    import numpy as np

    phi = np.deg2rad(strike)
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    m_rr = (np.sin(2.0 * delta) * np.sin(lambd)) * M0 / 1e14

    m_pp = (
        (
            np.sin(delta) * np.cos(lambd) * np.sin(2.0 * phi)
            - np.sin(2.0 * delta) * np.cos(phi) ** 2.0 * np.sin(lambd)
        )
        * M0
        / 1e14
    )

    m_tt = (
        (
            -np.sin(delta) * np.cos(lambd) * np.sin(2.0 * phi)
            - np.sin(2.0 * delta) * np.sin(phi) ** 2.0 * np.sin(lambd)
        )
        * M0
        / 1e14
    )

    m_rp = (
        (
            -np.cos(phi) * np.sin(lambd) * np.cos(2.0 * delta)
            + np.cos(delta) * np.cos(lambd) * np.sin(phi)
        )
        * M0
        / 1e14
    )

    m_rt = (
        (
            -np.sin(lambd) * np.sin(phi) * np.cos(2.0 * delta)
            - np.cos(delta) * np.cos(lambd) * np.cos(phi)
        )
        * M0
        / 1e14
    )

    m_tp = (
        (
            -np.sin(delta) * np.cos(lambd) * np.cos(2.0 * phi)
            - np.sin(2.0 * delta) * np.sin(2.0 * phi) * np.sin(lambd) / 2.0
        )
        * M0
        / 1e14
    )

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


def Direct(
    fwd: _Forward._AbstractForward,
    misfit: _Misfit._AbstractMisfit,
    event: obspy.core.event.Event,
    rec: instaseis.Receiver,
    phases: [str],
    components: [str],
    t_pre: [str],
    t_post: [str],
    depths: [float],
    phase_corrs: [float] = None,
    tstars: _Union[_List[float], _List[str]] = None,
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = False,
    list_to_correct_M0: [str] = None,
    output_folder: str = None,
    plot: bool = False,
    plot_extra_phases: [str] = None,
    color_plot: str = None,
    Ylims: [float] = None,
):
    print(f"Running direct inversion with model: {fwd.name}")
    print(f"and with {misfit.description}")
    M0 = 1e14

    if tstars is None:
        tstars = [None] * len(phases)

    if phase_corrs is None:
        phase_corrs = [0] * len(phases)

    if (fmin == None) or (fmax == None):
        print("Data will not be filtered due to fmin or fmax equal to None")
        filter_par = False
    else:
        filter_par = True

    if output_folder is None:
        output_folder = "."

    # TODO: IMPLEMENT LQT COORDINATE SYSTEM
    LQT_value = False
    baz = None
    inc = None

    """ PRE-PROCESS THE OBSERVED DATA """
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
    st_obs, sigmas = _PreProcess.prepare_event_data(
        event=event,
        phases=phases,
        components=components,
        slice=True,
        tts=obs_tt,
        t_pre=t_pre,
        t_post=t_post,
        filter=filter_par,
        fmin=fmin,
        fmax=fmax,
        zerophase=zerophase,
        noise_level=misfit.noise_level,
    )

    rec_in = instaseis.Receiver(
        latitude=90.0 - event.distance,
        longitude=0.0,
        network="XB",
        station="ELYSE",
        location="02",
    )

    for depth in depths:
        print(depth)
        ## Do inversion
        syn_tts = []
        for i, phase in enumerate(phases):
            syn_tt = fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance)
            syn_tts.append(syn_tt)

            syn_GF = fwd.get_greens_functions(
                comp=components[i],
                depth=depth,
                distance=event.distance,
                lat_src=90.0,
                lon_src=0.0,
                rec=rec_in,
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

            G = fwd.generate_G_matrix(
                st_GF=syn_GF,
                az=event.az,
                comp=components[i],
                slice=True,
                tt=syn_tt,
                t_pre=t_pre[i],
                t_post=t_post[i],
            )

            ## Weight matrix:
            start_weight = misfit.weights[i][0]
            end_weight = misfit.weights[i][1]

            samps = int(misfit.start_weight_len / misfit.dt)
            d_weight = _np.zeros_like(st_obs[i].data)
            d_weight[:samps] = start_weight
            d_weight[samps:] = end_weight

            Wd = 1 / (sigmas[i] ** 2 * d_weight)

            if i == 0:
                G_tot = G
                d_tot = st_obs[i].data
                Wd_tot = Wd

            else:
                G_tot = _np.vstack((G_tot, G))
                d_tot = _np.hstack((d_tot, st_obs[i].data))
                Wd_tot = _np.hstack((Wd_tot, Wd))

        Wd_tot = _np.diag(Wd_tot)

        # ---- Solve ----
        # TODO: Add weight matrix when solving the inverse
        A = G_tot.T @ Wd_tot @ G_tot
        B = G_tot.T @ Wd_tot @ d_tot

        try:
            M = _np.linalg.solve(A, B)
        except _np.linalg.LinAlgError:
            M = _np.linalg.lstsq(A, B)[0]

        # --- transform to r, theta, phi system ---
        mxx = M[0]
        myy = M[1]
        mzz = -mxx - myy
        mxy = M[2]
        mxz = M[3]
        myz = M[4]

        m_rr = mzz
        m_pp = myy
        m_tt = mxx
        m_rp = -myz
        m_rt = mxz
        m_tp = -mxy

        MT = [m_rr, m_pp, m_tt, m_rp, m_rt, m_tp]
        M0 = (
            m_rr ** 2 + m_tt ** 2 + m_pp ** 2 + 2 * m_rt ** 2 + 2 * m_rp ** 2 + 2 * m_tp ** 2
        ) ** 0.5 * 0.5 ** 0.5
        print("Scalar Moment: %.4e" % M0)
        MW = 2.0 / 3.0 * (_np.log10(M0) - 9.1)
        print("Magnitude: %.2f" % MW)

        ## Decompose MT into CLVD & DC:
        M_tensor = _np.array(
            [[MT[2], -MT[5], MT[4]], [-MT[5], MT[1], -MT[3]], [MT[4], -MT[3], MT[0]]]
        )
        M_CLVD, M_DC, F = _MTDecompose.Get_CLVD_DC(M_tensor)

        DC_MT = [M_DC[2, 2], M_DC[1, 1], M_DC[0, 0], -M_DC[1, 2], M_DC[0, 2], -M_DC[0, 1]]
        M0_DC = (
            M_DC[2, 2] ** 2
            + M_DC[0, 0] ** 2
            + M_DC[1, 1] ** 2
            + 2 * M_DC[0, 2] ** 2
            + 2 * M_DC[1, 2] ** 2
            + 2 * M_DC[0, 1] ** 2
        ) ** 0.5 * 0.5 ** 0.5

        CLVD_MT = [
            M_CLVD[2, 2],
            M_CLVD[1, 1],
            M_CLVD[0, 0],
            -M_CLVD[1, 2],
            M_CLVD[0, 2],
            -M_CLVD[0, 1],
        ]
        M0_CLVD = (
            M_CLVD[2, 2] ** 2
            + M_CLVD[0, 0] ** 2
            + M_CLVD[1, 1] ** 2
            + 2 * M_CLVD[0, 2] ** 2
            + 2 * M_CLVD[1, 2] ** 2
            + 2 * M_CLVD[0, 1] ** 2
        ) ** 0.5 * 0.5 ** 0.5

        ## Compute synthetics with the calculated MT
        syn_GFs = []
        st_syn = obspy.Stream()
        for i, phase in enumerate(phases):
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
            tr_syn = fwd.generate_synthetic_data(
                st_GF=syn_GF,
                focal_mech=DC_MT,
                M0=1.0,
                slice=True,
                tt=syn_tts[i],
                t_pre=t_pre[i],
                t_post=t_post[i],
            )
            st_syn += tr_syn
        ## Calculate the misfit
        chi = misfit.run_misfit(phases=phases, st_obs=st_obs, st_syn=st_syn, sigmas=sigmas)

        ## Calculate take-off angles P,S & pP
        takeoff_angles = ["P", "S", "pP"]
        angles = []
        for phase in takeoff_angles:
            angles.append(
                fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance, takeoffs=True)
            )

        """ Open .h5 file """
        file_name = (
            f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.hdf5"
        )
        f = _h5.File(pjoin(output_folder, file_name), "w")
        data_len = 5 + 3 * 6 + len(angles) + len(phases)
        file_len = 1
        f.create_dataset("samples", (file_len, data_len))
        """ Write into file """
        f["samples"][0, :] = [depth, F, M0, M0_DC, M0_CLVD] + MT + DC_MT + CLVD_MT + angles + chi

        if plot:
            # TODO: make this plot routine as a function
            MT = f["samples"][0, 5 : 5 + 6]
            DC_MT = f["samples"][0, 5 + 6 : 5 + 2 * 6]
            CLVD_MT = f["samples"][0, 5 + 6 : 5 + 2 * 6]

            angles = f["samples"][0, 5 + 3 * 6 : -len(phases)]
            angle_names = takeoff_angles
            chi = f["samples"][0, -len(phases) :]

            F = f["samples"][0, 1]
            M0 = f["samples"][0, 2]
            M0_DC = f["samples"][0, 3]
            M0_CLVD = f["samples"][0, 4]

            if color_plot is None:
                color_plot = "red"

            # For the plotting you need to re-arrange the order of the Moment tensor#
            # , since obspy uses a different order
            FULL = _np.array([MT[0], MT[2], MT[1], MT[4], MT[3], MT[5]])
            DC = _np.array([DC_MT[0], DC_MT[2], DC_MT[1], DC_MT[4], DC_MT[3], DC_MT[5]])
            CLVD = _np.array(
                [CLVD_MT[0], CLVD_MT[2], CLVD_MT[1], CLVD_MT[4], CLVD_MT[3], CLVD_MT[5]]
            )

            fig = _PostProcessing.Plot_Direct_BB(
                MT_Full=FULL / M0,
                Eps=F,
                MT_DC=DC / M0_DC,
                M0_DC=M0_DC,
                MT_CLVD=CLVD / M0_CLVD,
                M0_CLVD=M0_CLVD,
                azimuths=[event.az, event.az, event.az],
                inc_angles=angles,
                phase_names=angle_names,
                color=color_plot,
                height=19.0,
                horizontal=True,
            )

            plt.savefig(
                pjoin(
                    output_folder,
                    f"Direct_BB_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}.svg",
                ),
                dpi=300,
            )
            plt.close()

            MT = _np.expand_dims(DC_MT, axis=0)
            M0 = _np.expand_dims(M0_DC, axis=0)

            extra_arrs = []
            for j, extraphase in enumerate(plot_extra_phases):
                arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
                extra_arrs.append(arr)

            fig = _PostProcessing.waveform_plot(
                syn_GFs=syn_GFs,
                syn_tts=syn_tts,
                obs_tts=obs_tt,
                fwd=fwd,
                misfit_weight_len=misfit.start_weight_len,
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
            plt.savefig(
                pjoin(
                    output_folder,
                    f"Direct_waveforms_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}.svg",
                ),
                dpi=300,
            )
            plt.close()
        f.close()


def MH(self, event: obspy.core.event.Event):
    pass

