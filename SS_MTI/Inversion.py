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
    output_folder=None,
    plot=False,
):
    """
    Grid search over strike, dip, rake angles
    :param event: Obspy.event including waveforms and phase arrivals
    :param phases: list of phases to include in the inversion
    """
    print(f"Running grid search with model: {fwd.name}")
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

    for depth in depths:
        """ Open .h5 file """
        file_name = f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}.hdf5"
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

                    print(focal_mech)
                    print(_np.sum(chi))
                    print(M0_corr * M0)

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

            """ Extra phases to plot:"""
            extra_phases = ["PP", "SS", "pP", "sP", "PPP", "SSS"]

            sum_misfits = _np.sum(f["samples"][:, -len(phases) :], axis=1)
            nlowest = 3
            lowest_indices = sum_misfits.argsort()[0:nlowest]
            sdrs_total = f["samples"][:, 1:4]
            sdrs = sdrs_total[lowest_indices, :]
            M0_corrs_total = f["samples"][:, 5]
            M0_corrs = M0_corrs_total[lowest_indices]
            """ Beachball plot """

            fig = _PostProcessing.Plot_GS_BB(
                sdrs[:, 0],
                sdrs[:, 1],
                sdrs[:, 2],
                azimuths=[event.az, event.az, event.az],
                inc_angles=angles,
                phase_names=takeoff_angles,
                color="blue",
            )
            plt.savefig(
                pjoin(output_folder, f"GS_BBB_{depth}_{strike}_{dip}_{rake}_{misfit.name}.pdf")
            )
            plt.close()

            """ Waveform plot """
            fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(8, 6))

            for i, phase in enumerate(phases):
                for n in range(len(lowest_indices)):
                    tr_syn_full = fwd.generate_synthetic_data(
                        st_GF=syn_GFs[i], focal_mech=sdrs[n], M0=M0, slice=False,
                    )

                    tr_slice = tr_syn_full.slice(
                        starttime=fwd.or_time + syn_tts[i] - t_pre[i],
                        endtime=fwd.or_time + syn_tts[i] + t_post[i],
                    )

                    if n == 0:
                        st_obs_full, sigmas = _PreProcess.prepare_event_data(
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
                        ax[i].plot(
                            st_obs_full[i].times() - obs_tt[i], st_obs_full[i].data, lw=2, c="k",
                        )
                        ax[i].plot(
                            st_obs[i].times() - t_pre[i],
                            st_obs[i].data,
                            lw=4,
                            c="k",
                            label="Observed",
                        )
                        ax[i].plot(
                            tr_slice.times() - t_pre[i],
                            tr_slice.data * M0_corrs[n],
                            lw=4,
                            c="r",
                            label="Synthetic",
                        )
                    else:
                        ax[i].plot(
                            tr_slice.times() - t_pre[i], tr_slice.data * M0_corrs[n], lw=4, c="r"
                        )
                    ax[i].plot(
                        tr_syn_full.times() - (syn_tts[i] - fwd.start_cut),
                        tr_syn_full.data * M0_corrs[n],
                        lw=2,
                        c="r",
                    )
                    # ax[i].legend()

                    st = obspy.Stream()
                    st += tr_slice
                    st += st_obs[i]
                    if n == len(lowest_indices) - 1:
                        global_max = max([tr.data.max() for tr in st]) * 1.2
                        global_min = min([tr.data.min() for tr in st]) * 1.2
                        ax[i].set_ylim(global_min, global_max)
                        ax[i].axvline(x=t_post[i], c="grey", ls="dashed")
                        ax[i].axvline(x=-t_pre[i], c="grey", ls="dashed")
                        ax[i].axvspan(
                            -t_pre[i], misfit.start_weight_len, facecolor="grey", alpha=0.2
                        )
                        ax[i].axvline(x=0.0, c="grey")
                        ax[i].text(
                            0 + 0.1,
                            global_max * 0.8,
                            phase,
                            verticalalignment="center",
                            color="grey",
                            fontsize=6,
                        )
                        ax[i].text(
                            s="%s%s" % (phases[i], components[i]),
                            x=0.98,
                            y=0.7,
                            ha="right",
                            transform=ax[i].transAxes,
                            color="blue",
                            fontsize=20,
                        )

                        # Extra phase arrivals:
                        for j, extraphase in enumerate(extra_phases):
                            arr = fwd.get_phase_tt(
                                phase=extraphase, depth=depth, distance=event.distance
                            )
                            ax[i].axvline(x=arr - syn_tts[i], c="grey")
                            ax[i].text(
                                arr - syn_tts[i] + 0.1,
                                global_max * 0.8,
                                extraphase,
                                verticalalignment="center",
                                color="grey",
                                fontsize=6,
                            )
                        ax[i].get_yaxis().get_offset_text().set_visible(False)
                        ax_max = max(ax[i].get_yticks())
                        exponent_axis = _np.floor(_np.log10(ax_max)).astype(int)
                        ax[i].annotate(
                            r"$\times$10$^{%i}$" % (exponent_axis),
                            xy=(0.01, 0.82),
                            xycoords="axes fraction",
                        )

            fig.text(0.01, 0.5, "Displacement (m)", va="center", rotation="vertical", fontsize=18)
            fig.text(
                0.5,
                0.88,
                event.name,
                ha="center",
                va="bottom",
                size="x-large",
                color="blue",
                fontsize=18,
            )

            ax[0].legend(
                prop={"size": 12},
                loc="center left",
                bbox_to_anchor=(0.0, 0.95),
                bbox_transform=fig.transFigure,
            )

            ax[-1].set_xlim(-10.0, 60.0)
            ax[-1].set_xlabel("time after phase (s)", fontsize=18)
            plt.savefig(
                pjoin(
                    output_folder, f"GS_waveforms_{depth}_{strike}_{dip}_{rake}_{misfit.name}.pdf"
                )
            )
            plt.close()
        f.close()


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
    output_folder=None,
    plot=False,
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

        M = _np.linalg.solve(A, B)
        # M = _np.linalg.lstsq(A, B)[0]

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
        file_name = f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}.hdf5"
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
                color="red",
                height=19.0,
                horizontal=True,
            )

            plt.savefig(pjoin(output_folder, f"Direct_BB_{depth}_{misfit.name}.pdf"))
            plt.close()

            """ Extra phases to plot:"""
            extra_phases = ["PP", "SS", "pP", "sP", "PPP", "SSS"]

            fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(8, 6))
            for i, phase in enumerate(phases):
                """ Plot is generated with the DC solution """
                tr_syn_full = fwd.generate_synthetic_data(
                    st_GF=syn_GFs[i], focal_mech=DC_MT, M0=1.0, slice=False,
                )

                tr_slice = tr_syn_full.slice(
                    starttime=fwd.or_time + syn_tts[i] - t_pre[i],
                    endtime=fwd.or_time + syn_tts[i] + t_post[i],
                )

                st_obs_full, sigmas = _PreProcess.prepare_event_data(
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
                ax[i].plot(
                    st_obs_full[i].times() - obs_tt[i], st_obs_full[i].data, lw=2, c="k",
                )
                ax[i].plot(
                    st_obs[i].times() - t_pre[i], st_obs[i].data, lw=4, c="k", label="Observed",
                )
                ax[i].plot(
                    tr_slice.times() - t_pre[i], tr_slice.data, lw=4, c="r", label="Synthetic",
                )

                ax[i].plot(
                    tr_syn_full.times() - (syn_tts[i] - fwd.start_cut),
                    tr_syn_full.data,
                    lw=2,
                    c="r",
                )
                # ax[i].legend()

                st = obspy.Stream()
                st += tr_slice
                st += st_obs[i]
                global_max = max([tr.data.max() for tr in st]) * 1.2
                global_min = min([tr.data.min() for tr in st]) * 1.2
                ax[i].set_ylim(global_min, global_max)
                ax[i].axvline(x=0.0, c="grey")
                ax[i].axvline(x=t_post[i], c="grey", ls="dashed")
                ax[i].axvline(x=-t_pre[i], c="grey", ls="dashed")
                ax[i].text(
                    0 + 0.1,
                    global_max * 0.8,
                    phase,
                    verticalalignment="center",
                    color="grey",
                    fontsize=6,
                )

                ax[i].text(
                    s="%s%s" % (phases[i], components[i]),
                    x=0.98,
                    y=0.7,
                    ha="right",
                    transform=ax[i].transAxes,
                    color="blue",
                    fontsize=20,
                )

                # Extra phase arrivals:
                for j, extraphase in enumerate(extra_phases):
                    arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
                    ax[i].axvline(x=arr - syn_tts[i], c="grey")
                    ax[i].text(
                        arr - syn_tts[i] + 0.1,
                        global_max * 0.8,
                        extraphase,
                        verticalalignment="center",
                        color="grey",
                        fontsize=6,
                    )

                ax[i].get_yaxis().get_offset_text().set_visible(False)
                ax_max = max(ax[i].get_yticks())
                exponent_axis = _np.floor(_np.log10(ax_max)).astype(int)
                ax[i].annotate(
                    r"$\times$10$^{%i}$" % (exponent_axis),
                    xy=(0.01, 0.9),
                    xycoords="axes fraction",
                )
                ax[i].axvspan(-t_pre[i], misfit.start_weight_len, facecolor="grey", alpha=0.2)

            fig.text(0.01, 0.5, "Displacement (m)", va="center", rotation="vertical", fontsize=18)
            fig.text(
                0.5,
                0.88,
                event.name,
                ha="center",
                va="bottom",
                size="x-large",
                color="blue",
                fontsize=18,
            )

            ax[0].legend(
                prop={"size": 12},
                loc="center left",
                bbox_to_anchor=(0.0, 0.95),
                bbox_transform=fig.transFigure,
            )
            ax[-1].set_xlim(-10.0, 60.0)
            ax[-1].set_xlabel("time after phase (s)", fontsize=18)
            plt.savefig(pjoin(output_folder, f"Direct_waveforms_{depth}_{misfit.name}.pdf"))
            plt.close()
        f.close()


def MH(self, event: obspy.core.event.Event):
    pass

