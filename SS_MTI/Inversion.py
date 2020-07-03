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


def Grid_Search_run(
    fwd: _Forward._AbstractForward,
    misfit: _Misfit._AbstractMisfit,
    event: obspy.core.event.Event,
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
            syn_GF, syn_tt = fwd.get_greens_functions(
                phase=phase,
                comp=components[i],
                depth=depth,
                distance=event.distance,
                lat_src=event.latitude,
                lon_src=event.longitude,
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
                    f["samples"][iteration, :] = [depth, strike, dip, rake, M0, M0_corr] + chi
                    iteration += 1

                    print(focal_mech)
                    print(_np.sum(chi))
                    print(M0_corr * M0)

        if plot:
            sum_misfits = _np.sum(f["samples"][:, -len(phases) :], axis=1)
            nlowest = 1
            lowest_indices = sum_misfits.argsort()[0:nlowest]
            sdrs = f["samples"][lowest_indices, 1:4]
            M0_corrs = f["samples"][lowest_indices, 5]
            fig, ax = plt.subplots(nrows=len(phases), ncols=1, sharex="all", figsize=(8, 6))
            for n in range(len(lowest_indices)):
                for i, phase in enumerate(phases):
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
                    ax[i].legend()
                    ax[i].set_ylim(-5e-10, 5e-10)
                    ax[i].axvline(x=0.0, c="grey")
            ax[-1].set_xlim(-10.0, 60.0)
            plt.savefig(pjoin(output_folder, f"{depth}_{strike}_{dip}_{rake}.pdf"))
        f.close()


def Direct(self, event: obspy.core.event.Event):
    pass


def MH(self, event: obspy.core.event.Event):
    pass

