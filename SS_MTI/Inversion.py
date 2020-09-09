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
from obspy.signal.cross_correlation import xcorr_max, correlate


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
    sigmas = [i ** 2 for i in sigmas_noise]

    for depth in depths:
        M0_corrs_range = [1.26126e14]  # _np.linspace(0, 2, 1000)
        print(depth)
        """ Open .h5 file """
        file_name = f"GS_{event.name}_{depth}_{fmin}_{fmax}_{misfit.name}_{fwd.veloc_name}.hdf5"
        f = _h5.File(pjoin(output_folder, file_name), "w")
        data_len = 6 + len(phases)
        file_len = len(M0_corrs_range)  # len(strikes) * len(dips) * len(rakes)
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
                    for M0_corr_i in M0_corrs_range:

                        focal_mech = [strike, dip, rake]
                        st_syn = obspy.Stream()
                        st_syn_full = obspy.Stream()
                        misfit_amp = []

                        """ Generate the synthetic data"""
                        for i, phase in enumerate(phases):
                            # tr_syn = fwd.generate_synthetic_data(
                            #     st_GF=syn_GFs[i],
                            #     focal_mech=focal_mech,
                            #     M0=M0,
                            #     slice=True,
                            #     tt=syn_tts[i],
                            #     t_pre=t_pre[i],
                            #     t_post=t_post[i],
                            # )
                            tr_syn_full = fwd.generate_synthetic_data(
                                st_GF=syn_GFs[i], focal_mech=focal_mech, M0=M0, slice=False,
                            )

                            tr_syn = tr_syn_full.slice(
                                starttime=fwd.or_time + syn_tts[i] - t_pre[i],
                                endtime=fwd.or_time + syn_tts[i] + t_post[i],
                            )

                            if phases[i] + components[i] in list_to_correct_M0:
                                start_weight = misfit.weights[i][0]
                                end_weight = misfit.weights[i][1]

                                samps = int(misfit.start_weight_len / misfit.dt)
                                d_weight = _np.zeros_like(st_obs[i].data)
                                d_weight[:samps] = start_weight
                                d_weight[samps:] = end_weight

                                W = 1 / (d_weight)
                                # W = _np.diag(1 / (d_weight))

                                d_obs = _np.expand_dims(st_obs[i].data * W, axis=1)
                                d_syn = _np.expand_dims(tr_syn.data * W, axis=1)

                                amplitude = ((d_obs.T @ d_syn) / (d_obs.T @ d_obs))[0][0]
                                # print(amplitude)
                                misfit_amp.append(_np.abs(amplitude))
                                # misfit_amp.append(
                                #     (_np.sum(_np.abs(st_obs[i].data)))
                                #     / (_np.sum(_np.abs(tr_syn.data)))
                                # )
                                # misfit_amp.append(
                                #     (max(abs(st_obs[i].data)))
                                #     / (max(abs(tr_syn.data)))
                                # )

                            st_syn += tr_syn
                            st_syn_full += tr_syn_full

                        """ Multiply the data with the M0 correction"""
                        # M0_corr = 1e14 * M0_corr_i  #
                        M0_corr = M0_corr_i
                        # M0_corr = 1 / _np.mean(misfit_amp)
                        print(M0_corr)
                        # M0_corr = _np.sum(misfit_amp) / len(misfit_amp)  # 9.18202e12
                        # M0_corr = _np.exp(abs(_np.log(misfit_amp[0] /
                        #                                    misfit_amp[1])))

                        shift_CC = _np.zeros(len(phases))
                        misfit_CC = _np.zeros(len(phases))
                        shifts = {"P": None, "S": None}
                        for i, tr in enumerate(st_syn):
                            tr.data = tr.data * M0_corr

                            # st_syn_full[i].data = st_syn_full[i].data * M0_corr

                            # corrarray = correlate(tr.data, st_obs[i].data, domain="time", shift=40)
                            # shift_CC[i], misfit_CC[i] = xcorr_max(corrarray, abs_max=False)

                            # if shifts[phase] is None:
                            #     shifts[phase] = shift_CC[i]
                            # else:
                            #     # misfit_CC[iphase] = \
                            #     #     corrarray[(len(corrarray) - 1) // 2
                            #     #               + int(shifts[phases[iphase]])]
                            #     shift_CC[i] = shifts[phase]

                            # start_sample = int((syn_tts[i] - t_pre[i] - fwd.start_cut) / fwd.dt) + int(
                            #     shifts[phase]
                            # )
                            # end_sample = start_sample + len(tr.data)
                            # tr_shifted_syn = st_syn_full[i].data[start_sample:end_sample]

                            # dat = _np.vstack((tr.data,st_obs[i].data))
                            # with open(pjoin(output_folder, f"GS_{tr_syn.stats.channel}_{phases[i]}.txt"), 'wb') as file:
                            #     _np.save(file, dat, allow_pickle=False)
                            #     file.close()

                        """ Determine the misfit between syntetic and observed"""
                        chi = misfit.run_misfit(
                            phases=phases, st_obs=st_obs, st_syn=st_syn, sigmas=sigmas
                        )
                        # print(chi)
                        """ Write into file"""
                        f["samples"][iteration, :] = [depth, strike, dip, rake, M0, M0_corr] + chi
                        iteration += 1

                        # print(focal_mech)
                        # print(_np.sum(chi))
                        # print(M0_corr * M0)

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
            lowest_indices = lowest_ind[inds]
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
                pjoin(
                    output_folder,
                    f"GS_waveforms_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}.svg",
                ),
                dpi=300,
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
    output_folder: str = None,
    plot: bool = False,
    plot_extra_phases: [str] = None,
    color_plot: str = None,
    Ylims: [float] = None,
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
    sigmas = [i ** 2 for i in sigmas_noise]
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

            Wd = 1 / (sigmas[i] * d_weight)
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
            # with open(pjoin(output_folder, f"Direct_{tr_syn.stats.channel}_{phase}.txt"), 'wb') as file:
            #     _np.save(file,dat, allow_pickle=False)
            #     file.close()

        ## Calculate the misfit
        chi = misfit.run_misfit(phases=phases, st_obs=st_obs, st_syn=st_syn, sigmas=sigmas)
        # print(chi)
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
            CLVD_MT = f["samples"][0, 5 + 2 * 6 : 5 + 3 * 6]

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
                pjoin(
                    output_folder,
                    f"Direct_BB_{event.name}_{depth}_{misfit.name}_{fwd.veloc_name}.svg",
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

