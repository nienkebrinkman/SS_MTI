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
from os.path import exists, isfile, join
from os import listdir, makedirs
from obspy.taup import TauPyModel
from obspy import UTCDateTime as utct
from numpy import linalg as _LA
from typing import List as _List, Union as _Union
from obspy.signal.cross_correlation import xcorr_max, correlate
from scipy.optimize import approx_fprime as _af
import mpi4py.MPI


from SS_MTI import PreProcess as _PreProcess
from SS_MTI import Forward as _Forward
from SS_MTI import Misfit as _Misfit
from SS_MTI import MTDecompose as _MTDecompose
from SS_MTI import PostProcessing as _PostProcessing
from SS_MTI import Gradient as _Gradient


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
    Parallel: bool = False,
):
    """
    Grid search over strike, dip, rake angles
    :param event: Obspy.event including waveforms and phase arrivals
    :param phases: list of phases to include in the inversion
    """
    print(f"Running grid search with model: {fwd.name}")
    print(f"and with {misfit.description}")
    M0 = 1

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
    st_obs, sigmas_noise = _PreProcess.prepare_event_data(
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

    # # TODO: remove this again:
    # if event.name == "Test_Event":
    #     sigmas_noise[0] = 1.7711953652440284e-11
    #     sigmas_noise[3] = 4.5996573530998017e-11

    # sigmas_noise[2] = sigmas_noise[1]
    # sigmas_noise[3] = sigmas_noise[0]
    # sigmas_noise[4] = sigmas_noise[1]
    # print(sigmas)
    # sigmas_model = [9e-10, 2e-9, 1e-9, 7e-10, 2e-9]

    # sigmas = [i ** 2 + j ** 2 for i, j in zip(sigmas_noise, sigmas_model)]
    variances = [i ** 2 for i in sigmas_noise]

    if Parallel:
        rank = mpi4py.MPI.COMM_WORLD.Get_rank()
        size = mpi4py.MPI.COMM_WORLD.Get_size()

    for iPar, depth in enumerate(depths):
        if Parallel:
            if iPar % size != rank:
                continue
            print(f"Depth number {depth} being done by processor {rank} of {size}")
        # M0_corrs_range = [1.26126e14]  # _np.linspace(0, 2, 1000)
        print(depth)
        """ Open .h5 file """
        file_name = f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}_{event.baz}.hdf5"
        f = _h5.File(join(output_folder, file_name), "w")
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

        shift_CC = _np.zeros(len(phases))
        misfit_CC = _np.zeros(len(phases))

        for strike in strikes:
            for dip in dips:
                for rake in rakes:
                    shifts = {"P": None, "S": None}

                    focal_mech = [strike, dip, rake]
                    st_syn = obspy.Stream()
                    st_syn_full = obspy.Stream()
                    misfit_amp = []

                    """ Generate the synthetic data"""
                    for i, phase in enumerate(phases):
                        tr_syn_full = fwd.generate_synthetic_data(
                            st_GF=syn_GFs[i], focal_mech=focal_mech, M0=M0, slice=False,
                        )

                        tr_syn = tr_syn_full.slice(
                            starttime=fwd.or_time + syn_tts[i] - t_pre[i],
                            endtime=fwd.or_time + syn_tts[i] + t_post[i],
                        )

                        # """ Calculate cross-correlation """
                        # from obspy.signal.cross_correlation import xcorr_max, correlate

                        # max_shift = int(1.5 / fwd.dt)
                        # corrarray = correlate(tr_syn, st_obs[i], domain="time", shift=max_shift)
                        # shift_CC[i], misfit_CC[i] = xcorr_max(corrarray, abs_max=False)

                        # if shifts[phases[i]] is None:
                        #     shifts[phases[i]] = shift_CC[i]
                        #     # print(shift_CC[iphase], misfit_CC[iphase],
                        #     #       corrarray[(len(corrarray) - 1) // 2
                        #     #                 + int(shifts[phases[iphase]])])
                        # else:
                        #     misfit_CC[i] = corrarray[
                        #         (len(corrarray) - 1) // 2 + int(shifts[phases[i]])
                        #     ]
                        #     shift_CC[i] = shifts[phases[i]]

                        # shift_in_sec = shift_CC[i] * fwd.dt

                        # tr_syn = tr_syn_full.slice(
                        #     starttime=fwd.or_time + syn_tts[i] - t_pre[i] - shift_in_sec,
                        #     endtime=fwd.or_time + syn_tts[i] + t_post[i] - shift_in_sec,
                        # )

                        if phases[i] + components[i] in list_to_correct_M0:
                            d_obs = _np.expand_dims(st_obs[i].data, axis=1)
                            d_syn = _np.expand_dims(tr_syn.data, axis=1)

                            ## Weight matrix:
                            start_weight = misfit.weights[i][0]
                            end_weight = misfit.weights[i][1]

                            samps = int(misfit.start_weight_len / misfit.dt)
                            d_weight = _np.zeros_like(st_obs[i].data)
                            d_weight[:samps] = start_weight
                            d_weight[samps:] = end_weight

                            Wd = 1 / (variances[i] * d_weight)
                            # Wd = 1 / (_np.std(st_obs[i].data) ** 2 * d_weight)

                            if i == 0:
                                G_tot = d_syn
                                d_tot = d_obs
                                Wd_tot = Wd

                            else:
                                G_tot = _np.vstack((G_tot, d_syn))
                                d_tot = _np.vstack((d_tot, d_obs))
                                Wd_tot = _np.hstack((Wd_tot, Wd))

                        st_syn += tr_syn
                        st_syn_full += tr_syn_full

                    Wd_tot = _np.diag(Wd_tot)

                    # ---- Solve ----
                    A = G_tot.T @ Wd_tot @ G_tot
                    B = G_tot.T @ Wd_tot @ d_tot

                    M0_corr = _np.abs(_np.linalg.solve(A, B)[0][0])

                    """ Multiply the data with the M0 correction"""
                    shift_CC = _np.zeros(len(phases))
                    misfit_CC = _np.zeros(len(phases))
                    shifts = {"P": None, "S": None}
                    for i, tr in enumerate(st_syn):
                        tr.data = tr.data * M0_corr

                    """ Determine the misfit between syntetic and observed"""
                    chi = misfit.run_misfit(
                        phases=phases, st_obs=st_obs, st_syn=st_syn, variances=variances
                    )
                    # print(chi)
                    """ Write into file"""
                    f["samples"][iteration, :] = [depth, strike, dip, rake, M0, M0_corr] + chi
                    iteration += 1

        if plot:
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

            Total_L2_GS = _np.sum(f["samples"][:, -len(phases) :], axis=1)
            lowest_ind = Total_L2_GS.argsort()
            Total_L2_GS.sort()
            misfit_low = Total_L2_GS[:] - Total_L2_GS[0]
            uncert = 0.05 * Total_L2_GS[0]
            inds = _np.where(misfit_low < uncert)
            lowest_indices = lowest_ind[inds][:10]
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
                join(
                    output_folder,
                    f"GS_BBB__{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}_{event.baz}.svg",
                ),
                dpi=300,
            )
            plt.close()

            if plot_extra_phases is not None:
                extra_arrs = []
                for j, extraphase in enumerate(plot_extra_phases):
                    arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
                    extra_arrs.append(arr)
            else:
                extra_arrs = None

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

            plt.savefig(
                join(
                    output_folder,
                    f"GS_waveforms_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}_{event.baz}.svg",
                ),
                dpi=300,
            )
            plt.close()
        f.close()
    if Parallel:
        mpi4py.MPI.COMM_WORLD.Barrier()


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
    Parallel: bool = False,
):
    print(f"Running direct inversion with model: {fwd.name}")
    print(f"and with {misfit.description}")

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
    st_obs, sigmas_noise = _PreProcess.prepare_event_data(
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

    # # TODO: remove this again:
    # if event.name == "Test_Event":
    #     sigmas_noise[0] = 1.7711953652440284e-11
    #     sigmas_noise[3] = 4.5996573530998017e-11
    # print(sigmas)
    # sigmas[2] = sigmas[1]
    # sigmas[3] = sigmas[0]
    # sigmas[4] = sigmas[1]

    # sigmas_model = [9e-10, 2e-9, 1e-9, 7e-10, 2e-9]

    # sigmas = [i ** 2 + j ** 2 for i, j in zip(sigmas_noise, sigmas_model)]
    variances = [i ** 2 for i in sigmas_noise]
    rec_in = instaseis.Receiver(
        latitude=90.0 - event.distance,
        longitude=0.0,
        network="XB",
        station="ELYSE",
        location="02",
    )

    if Parallel:
        rank = mpi4py.MPI.COMM_WORLD.Get_rank()
        size = mpi4py.MPI.COMM_WORLD.Get_size()

    for iPar, depth in enumerate(depths):
        if Parallel:
            if iPar % size != rank:
                continue
            print(f"Depth number {depth} being done by processor {rank} of {size}")
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

            Wd = 1 / (variances[i] * d_weight)
            # Wd = 1 / (_np.std(st_obs[i].data) ** 2 * d_weight)

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
            print("least-square")
            M = _np.linalg.lstsq(A, B)[0]

        # --- Test condition number ---
        cond_nr = _LA.cond(G_tot.T @ Wd_tot @ G_tot)
        eigenValues, eigenVectors = _LA.eig(G_tot.T @ Wd_tot @ G_tot)

        idx = eigenValues.argsort()[::-1]
        Eigval = eigenValues[idx]
        Eigvec = eigenVectors[:, idx]

        max_eigval = Eigval.max()
        max_eigvec = Eigvec[:, _np.argmax(Eigval)]

        moment_names = ["mxx", "myy", "mxy", "mxz", "myz"]
        colors = ["k", "r", "green", "magenta", "darkorange"]
        fig, ax = plt.subplots(6, 1, sharex=True, figsize=(8, 10))
        ax[0].scatter(_np.arange(len(Eigval)), Eigval, color=colors, s=100)
        ax[0].set_yscale("log")
        ax[0].set_ylabel("Eigenvalues", fontsize=14)
        ax[0].tick_params(axis="both", which="major", labelsize=14)
        ax[0].tick_params(axis="both", which="minor", labelsize=10)
        ax[0].set_title("Condition number: {:.2f}".format(cond_nr), fontsize=20)
        for ax_i in range(len(Eigval)):
            ax_x = _np.arange(len(Eigvec[:, ax_i]))

            for ax_ii in range(len(ax_x)):
                if Eigvec[ax_ii, ax_i] > 0:
                    if ax_ii == 0:
                        ax[ax_i + 1].plot(
                            ax_x[ax_ii],
                            Eigvec[ax_ii, ax_i],
                            "k^",
                            label=f"Eigenvalue {ax_i}",
                            color=colors[ax_i],
                            markersize=10,
                        )
                    else:
                        ax[ax_i + 1].plot(
                            ax_x[ax_ii],
                            Eigvec[ax_ii, ax_i],
                            "k^",
                            color=colors[ax_i],
                            markersize=10,
                        )
                    ax[ax_i + 1].plot(
                        [ax_x[ax_ii], ax_x[ax_ii]], [0.0, Eigvec[ax_ii, ax_i]], c=colors[ax_i]
                    )
                else:
                    if ax_ii == 0:
                        ax[ax_i + 1].plot(
                            ax_x[ax_ii],
                            Eigvec[ax_ii, ax_i],
                            "kv",
                            label=f"Eigenvalue {ax_i}",
                            color=colors[ax_i],
                            markersize=10,
                        )
                    else:
                        ax[ax_i + 1].plot(
                            ax_x[ax_ii],
                            Eigvec[ax_ii, ax_i],
                            "kv",
                            color=colors[ax_i],
                            markersize=10,
                        )
                    ax[ax_i + 1].plot(
                        [ax_x[ax_ii], ax_x[ax_ii]], [Eigvec[ax_ii, ax_i], 0.0], c=colors[ax_i]
                    )

            ax[ax_i + 1].set_ylabel("Eigvector", fontsize=14)
            ax[ax_i + 1].legend(loc="upper right")
            ax[ax_i + 1].tick_params(axis="both", which="major", labelsize=14)
            ax[ax_i + 1].tick_params(axis="both", which="minor", labelsize=10)
            ax[ax_i + 1].set_ylim(-1.0, 1.0)

        ax[-1].xaxis.set_ticks(_np.arange(len(Eigval)))
        ax[-1].set_xticklabels(moment_names)
        ax[-1].set_xlabel("Parameters", fontsize=20)
        plt.tight_layout()
        plt.savefig(
            join(
                output_folder,
                f"Condition_nr_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}_{event.baz}.svg",
            ),
            dpi=600,
        )
        plt.close()

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

        # m_rr = -M[0] - M[1]
        # m_pp = M[1]
        # m_tt = M[0]
        # m_rp = M[4]
        # m_rt = M[3]
        # m_tp = M[2]

        MT = [m_rr, m_tt, m_pp, m_rt, m_rp, m_tp]
        M0 = (
            m_rr ** 2 + m_tt ** 2 + m_pp ** 2 + 2 * m_rt ** 2 + 2 * m_rp ** 2 + 2 * m_tp ** 2
        ) ** 0.5 * 0.5 ** 0.5  # TODO HAAKJES hallo
        print("Full Scalar Moment: %.4e" % M0)
        MW = 2.0 / 3.0 * (_np.log10(M0) - 9.1)
        print("Full magnitude: %.2f" % MW)

        ## Decompose MT into CLVD & DC:
        # USE = _np.array(
        #     [[MT[0], MT[4], MT[3]], [MT[4], MT[2], MT[5]], [MT[3], MT[5], MT[1]]]
        # )
        # USE_system  = _np.matrix([[0., -1., 0.], [0., 0., 1.], [-1., 0., 0.]],dtype=_np.float)
        # NED_system  = _np.matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        #                  dtype=_np.float)
        # r = USE_system.I * NED_system.I
        # NED = _np.dot(r, _np.dot(USE, r.I))
        # M_CLVD_NED, M_DC_NED, F = _MTDecompose.Get_CLVD_DC(NED)

        # r = NED_system.I * USE_system.I
        # M_CLVD_USE = _np.dot(r, _np.dot(M_CLVD_NED, r.I))
        # CLVD_MT = [M_CLVD_USE[0,0],M_CLVD_USE[2,2],M_CLVD_USE[1,1],M_CLVD_USE[0,2],M_CLVD_USE[0,1],M_CLVD_USE[1,2]]
        # M_DC_USE = _np.dot(r, _np.dot(M_DC_NED, r.I))
        # DC_MT = [M_DC_USE[0,0],M_DC_USE[2,2],M_DC_USE[1,1],M_DC_USE[0,2],M_DC_USE[0,1],M_DC_USE[1,2]]

        # M_tensor = _np.array(
        #     [[MT[2], -MT[5], MT[4]], [-MT[5], MT[1], -MT[3]], [MT[4], -MT[3], MT[0]]]
        # )
        M_tensor = _np.array(
            [[MT[1], -MT[5], MT[3]], [-MT[5], MT[2], -MT[4]], [MT[3], -MT[4], MT[0]]]
        )
        M_CLVD, M_DC, F = _MTDecompose.Get_CLVD_DC(M_tensor)

        DC_MT = [M_DC[2, 2], M_DC[0, 0], M_DC[1, 1], M_DC[0, 2], -M_DC[1, 2], -M_DC[0, 1]]
        M0_DC = (
            M_DC[2, 2] ** 2
            + M_DC[0, 0] ** 2
            + M_DC[1, 1] ** 2
            + 2 * M_DC[0, 2] ** 2
            + 2 * M_DC[1, 2] ** 2
            + 2 * M_DC[0, 1] ** 2
        ) ** 0.5 * 0.5 ** 0.5
        print("DC Scalar Moment: %.4e" % M0_DC)
        print("Epsilon value: %.2f" % F)

        CLVD_MT = [
            M_CLVD[2, 2],
            M_CLVD[0, 0],
            M_CLVD[1, 1],
            M_CLVD[0, 2],
            -M_CLVD[1, 2],
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
            # true=[0.0 ,-4.86706276927e+13 ,4.86706276927e+13 ,-0.0017206287528, 0.00298021642082 ,2.81e+13]
            tr_syn = fwd.generate_synthetic_data(
                st_GF=syn_GF,
                focal_mech=DC_MT,
                M0=1.0,
                slice=True,
                tt=syn_tts[i],
                t_pre=t_pre[i],
                t_post=t_post[i],
            )
            # tr_syn.data = _np.zeros_like(tr_syn.data)
            st_syn += tr_syn

            # moment_new = _np.expand_dims(
            #     [DC_MT[1], DC_MT[2], -DC_MT[5], DC_MT[3], -DC_MT[4]], axis=1
            # )
            # if i == 0:
            #     G_new = G_tot[: len(st_obs[i].data), :]
            # else:
            #     G_new = G_tot[
            #         len(st_obs[i - 1].data) : len(st_obs[i - 1].data) + len(st_obs[i].data), :
            #     ]
            # d_new = G_new @ moment_new

            # residual = tr_syn.data - d_new

            # dat = _np.vstack((tr_syn.data,st_obs[i].data))
            # with open(join(output_folder, f"Direct_{tr_syn.stats.channel}_{phase}.txt"), 'wb') as file:
            #     _np.save(file,dat, allow_pickle=False)
            #     file.close()

        ## Calculate the misfit
        chi = misfit.run_misfit(phases=phases, st_obs=st_obs, st_syn=st_syn, variances=variances)
        # print(chi)
        ## Calculate take-off angles P,S & pP
        takeoff_angles = ["P", "S", "pP"]
        angles = []
        for phase in takeoff_angles:
            angles.append(
                fwd.get_phase_tt(phase=phase, depth=depth, distance=event.distance, takeoffs=True)
            )

        """ Open .h5 file """
        file_name = f"Direct_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}_{event.baz}.hdf5"
        f = _h5.File(join(output_folder, file_name), "w")
        data_len = 6 + 3 * 6 + len(angles) + len(phases)
        file_len = 1
        f.create_dataset("samples", (file_len, data_len))
        """ Write into file """
        f["samples"][0, :] = (
            [depth, cond_nr, F, M0, M0_DC, M0_CLVD] + MT + DC_MT + CLVD_MT + angles + chi
        )

        if plot:
            # TODO: make this plot routine as a function
            MT = f["samples"][0, 6 : 6 + 6]
            DC_MT = f["samples"][0, 6 + 6 : 6 + 2 * 6]
            CLVD_MT = f["samples"][0, 6 + 2 * 6 : 6 + 3 * 6]

            angles = f["samples"][0, 6 + 3 * 6 : -len(phases)]
            angle_names = takeoff_angles
            chi = f["samples"][0, -len(phases) :]

            cond_nr = f["samples"][0, 1]
            F = f["samples"][0, 2]
            M0 = f["samples"][0, 3]
            M0_DC = f["samples"][0, 4]
            M0_CLVD = f["samples"][0, 5]

            if color_plot is None:
                color_plot = "red"

            # For the plotting you need to re-arrange the order of the Moment tensor#
            # , since obspy uses a different order
            FULL = MT  # _np.array([MT[0], MT[2], MT[1], MT[4], MT[3], MT[5]])
            DC = DC_MT  # _np.array([DC_MT[0], DC_MT[2], DC_MT[1], DC_MT[4], DC_MT[3], DC_MT[5]])
            CLVD = CLVD_MT  # _np.array(
            #     [CLVD_MT[0], CLVD_MT[2], CLVD_MT[1], CLVD_MT[4], CLVD_MT[3], CLVD_MT[5]]
            # )
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
                join(
                    output_folder,
                    f"Direct_BB_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}_{event.baz}.svg",
                ),
                dpi=300,
            )
            plt.close()

            MT = _np.expand_dims(DC_MT, axis=0)
            M0 = _np.expand_dims(M0_DC, axis=0)

            if plot_extra_phases is not None:
                extra_arrs = []
                for j, extraphase in enumerate(plot_extra_phases):
                    arr = fwd.get_phase_tt(phase=extraphase, depth=depth, distance=event.distance)
                    extra_arrs.append(arr)
            else:
                extra_arrs = None

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
                join(
                    output_folder,
                    f"Direct_waveforms_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}_{event.baz}.svg",
                ),
                dpi=300,
            )
            plt.close()
        f.close()
    if Parallel:
        mpi4py.MPI.COMM_WORLD.Barrier()


def gradient_descent(
    bin_path: str,
    save_path: str,
    epsilon: float,
    update_nr: int,
    dt: float,
    sigmas: [float],
    st_obs_w: obspy.Stream,
    current_update: int = 0,
    prior_crfl_filepath: str = None,
    alphas: [float] = [1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1],
    fmin: float = None,
    fmax: float = None,
    phases: [str] = ["P", "S", "P", "S", "S"],
    comps: [str] = ["Z", "T", "R", "Z", "R"],
    t_pres: [int] = [1, 1, 1, 1, 1],
    t_posts: [int] = [30, 30, 30, 30, 30],
):
    """ 
    This function will do a gradient-descent step based on a line-search
    :param bin_path: filepath to reflectivity binary
    :param save_path: path where all updates will be saved
    :param epsilon: step size w.r.t. each parameter
    :param update_nr: amount of updates you want to do
    :param dt: sampling rate [second/sample]
    :param sigmas: expected standard deviation for each phase
    :param st_obs_w: stream with windowed observed data
    :param current_update: current update nr, it restarts from previous update
    :param prior_crfl_filepath: only necessary when update_nr = 0
    :param alphas: gradient step sizes to test in the line search
    :param fmin: highpass frequency band
    :param fmax: lowpass frequency band
    :param phases: phases to window
    :param comps: components to window
    :param t_pres: length before arrival
    :param t_posts: length after arrival
    """

    #     assert (
    #         current_update == 0 and prior_crfl_filepath is not None
    #     ), "if current_update = 0, you have to specify filepath of your prior crfl.dat"

    save_path_OG = save_path

    if current_update != 0:
        """ Check where the previous update ended and take this crfl.dat file as prior file"""
        prev_update = current_update - 1
        prev_it = max(
            [
                int(f.strip("It_"))
                for f in listdir(join(save_path_OG, f"Update_{prev_update}"))
                if f.startswith("It_")
            ]
        )
        prior_crfl_filepath = join(
            save_path_OG, f"Update_{prev_update}", f"It_{prev_it}", "crfl.dat"
        )
        prev_m0 = [
            f
            for f in listdir(join(save_path, f"Update_{prev_update}"))
            if f.startswith("m1_")
            if isfile(join(save_path, f"Update_{prev_update}", f))
        ][0]
        m0 = _np.load(join(save_path_OG, f"Update_{prev_update}", prev_m0,))
    else:
        # NOTE: we have to convert fm parameters:
        # reflectivity fm: mtt,mtp,mrt,mpp,mrp,mrr
        # this code(also instaseis) fm: mrr,mtt,mpp,mrt,mrp,mtp
        with open(prior_crfl_filepath, "r") as f:
            data = f.readlines()
            fm = _np.array(data[-8].split(), dtype=float)
        m0 = _np.array(
            [fm[5], fm[0], fm[3], fm[2], -fm[4] + 0, -fm[1] + 0, float(data[9].split()[0])]
        )

    while current_update < update_nr:
        if not exists(join(save_path_OG, f"Update_{current_update}")):
            makedirs(join(save_path_OG, f"Update_{current_update}"))
        save_path = join(save_path_OG, f"Update_{current_update}")

        """ Calculating the gradient with given epsilon """
        src_str = _Gradient.SRC_STR(
            binary_file_path=bin_path,
            prior_dat_filepath=prior_crfl_filepath,
            save_folder=save_path,
            phases=phases,
            components=comps,
            t_pres=t_pres,
            t_posts=t_posts,
            depth=True,
            vpvs=False,
            fmin=fmin,
            fmax=fmax,
            dt=dt,
            sigmas=sigmas,
            zerophase=False,
            start_it=0,
        )

        dxi_dms = _np.zeros((len(m0), 1))
        if not isfile(join(save_path, "dxi_dms.npy")):
            dxi_dm = _af(
                m0,
                src_str.misfit,
                epsilon
                * _np.array(
                    [
                        _np.mean(m0[:-1]),
                        _np.mean(m0[:-1]),
                        _np.mean(m0[:-1]),
                        _np.mean(m0[:-1]),
                        _np.mean(m0[:-1]),
                        _np.mean(m0[:-1]),
                        0.1 * m0[-1],
                    ]
                ),
                st_obs_w,
            )
            dxi_dms[:, 0] = dxi_dm
            _np.save(join(save_path, "dxi_dms.npy"), dxi_dms)
        else:
            print("dxi_dms.npy already exists in this folder, reads in the existing file")
            dxi_dms = _np.load(join(save_path, "dxi_dms.npy"))

        """ Doing update using a line-search (i.e., update based on best alpha) """
        if not isfile(join(save_path, f"m1s_{epsilon}.npy")):
            m1s = _np.zeros((len(m0), len(alphas)))
            X1s = _np.zeros(len(alphas))
            for i, alpha in enumerate(alphas):
                m1s[:, i] = m0 - dxi_dm * alpha

                X1s[i] = src_str.misfit(m1s[:, i], st_obs_w)
            _np.save(join(save_path, f"m1s_{epsilon}.npy"), m1s)
            _np.save(join(save_path, f"X1s_{epsilon}.npy"), X1s)
        else:
            m1s = _np.load(join(save_path, f"m1s_{epsilon}.npy"))
            X1s = _np.load(join(save_path, f"X1s_{epsilon}.npy"))

        min_misfit = X1s.argmin()
        min_alpha = alphas[min_misfit]
        m1 = m1s[:, min_misfit]
        _np.save(join(save_path, f"misfit.npy"), X1s.min())
        _np.save(join(save_path, f"alpha.npy"), min_alpha)
        _np.save(join(save_path, f"m1_eps_{epsilon}_alpha_{min_alpha}.npy"), m1)

        """ Do a final forward run with the achieved m1 model: """
        st_m1 = src_str.forward(m1)
        st_m1.write(join(save_path, "st_m1.mseed"), format="MSEED")

        """
        Make the latest iteration in this update, which is based on m1,
        the new prior for the next update
        """
        update_it = src_str.it - 1
        print(f"this is the iteration used for next update: {update_it}")
        prior_crfl_filepath = join(
            save_path_OG, f"Update_{current_update}", f"It_{update_it}", "crfl.dat"
        )

        current_update += 1
        m0 = m1


def gauss_newton():

    J = 1

    J_inv = J.T @ J
    J_d = J @ d_tot

    try:
        M_upd = _np.linalg.solve(J_inv, J_d)
    except _np.linalg.LinAlgError:
        print("least-square")
        M_upd = _np.linalg.lstsq(J_inv, J_d)[0]
    m1s[:, i] = m0 - M_upd


def MH(self, event: obspy.core.event.Event):
    pass

