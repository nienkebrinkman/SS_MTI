#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Misfit method

:copyright:
    Nienke Brinkman (nienke.brinkman@erdw.ethz.ch), 2020
:license:
    None
"""


import abc
import obspy as _obspy
import numpy as _np
from typing import List as _List

# from obspy.signal.cross_correlation import xcorr_max as _xcorr_max, correlate as _correlate

from SS_MTI import PhaseTracer as _PhaseTracer
from SS_MTI import GreensFunctions as _GreensFunctions


class _AbstractMisfit(metaclass=abc.ABCMeta):
    name = "abstract Misfit"

    def __init__(self):
        pass

    @abc.abstractmethod
    def run_misfit(self):
        pass


class L2(_AbstractMisfit):
    name = "L2"
    description = "L2 based misfit"
    noise_level = True  # Do you need to determine level of noise for this misfit class?

    def __init__(self, weights: [_List[int]], start_weight_len: float, dt: float):
        """
        L2 misfit initializer
        :param weights: a list for each phase/comp combination that each have another list [start_weight, end_weight]
        :param start_weight_len: the amount of seconds that the start_weight is used in the window

        """
        self.weights = weights
        self.start_weight_len = start_weight_len
        self.dt = dt

    def run_misfit(
        self, phases: [str], st_obs: _obspy.Stream, st_syn: _obspy.Stream, sigmas: [float],
    ) -> [float]:
        """
        L2 misfit is calculated with a weighting matrix based on a noise array before arrival

        """
        assert (len(st_obs) == len(st_syn)) and (len(st_obs[0].data) == len(st_syn[0].data)), (
            "st_obs and st_syn should have equal amount"
            " of traces AND the trace lengths should be the same"
        )

        samps = int(self.start_weight_len / self.dt)

        assert (
            len(st_obs[0].data) > samps
        ), "The start_weight_len should not be longer than the window you are inverting for"

        misfit_L2 = []
        for i in range(len(phases)):
            d_obs = _np.expand_dims(st_obs[i].data, axis=1)
            d_syn = _np.expand_dims(st_syn[i].data, axis=1)

            start_weight = self.weights[i][0]
            end_weight = self.weights[i][1]

            d_weight = _np.zeros_like(st_obs[i].data)
            d_weight[:samps] = start_weight
            d_weight[samps:] = end_weight

            misfit_L2.append(
                0.5
                * (
                    (d_obs - d_syn).T
                    @ (_np.expand_dims(1 / (sigmas[i] ** 2 * d_weight), axis=1) * (d_obs - d_syn))
                )[0][0]
            )
        return misfit_L2


class CC(_AbstractMisfit):
    name = "CC"
    description = "cross-correlation based misfit"
    noise_level = False  # Do you need to determine level of noise for this misfit class?

    def __init__(self, shift_samples: float = 20):
        self.shift = shift_samples

    def run_misfit(
        self, phases: [str], st_obs: _obspy.Stream, st_syn: _obspy.Stream, sigmas: [float]
    ):
        assert (len(st_obs) == len(st_syn)) and (len(st_obs[0].data) == len(st_syn[0].data)), (
            "st_obs and st_syn should have equal amount"
            " of traces AND the trace lengths should be the same"
        )
        # TODO: FIX THE CROSS-CORRELATION FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        shift_CC = [None] * len(phases)
        misfit_CC = [None] * len(phases)

        for iphase in range(len(phases)):
            import obspy.signal.cross_correlation as cc
            import numpy as np

            """ nothing worked so I used a padd zero function, still cross-correlation doesnt work, no clue why"""
            offset = [25]
            a = st_obs[iphase].data
            resulta = np.zeros(len(st_obs[iphase].data) + 50)
            insertHere = [slice(offset[dim], offset[dim] + a.shape[dim]) for dim in range(a.ndim)]
            resulta[insertHere] = a

            b = st_syn[iphase].data
            resultb = np.zeros(len(st_obs[iphase].data) + 50)
            insertHere = [slice(offset[dim], offset[dim] + b.shape[dim]) for dim in range(a.ndim)]
            resultb[insertHere] = b

            cc_obspy = cc.correlate(resultb, resulta, len(resulta))

            max_shift = 1.0
            dt = st_obs[iphase].stats.delta
            max_shift_sample = int(max_shift / dt)

            cc_obspy = cc.correlate(st_syn[iphase].data, st_obs[iphase].data, max_shift_sample)
            shift_centered, MAX_CC = cc.xcorr_max(cc_obspy, abs_max=False)

            c = cc.correlate(
                st_syn[iphase].data, st_obs[iphase].data, len(st_syn[iphase].data)
            )  # a0 is the trace you want to shift in the end!!
            shift, Max_c = cc.xcorr_max(c, abs_max=False)

            corr_array = _correlate(
                st_syn[iphase].data, st_obs[iphase].data, domain="time", shift=self.shift,
            )
            shift_CC[iphase], misfit_CC[iphase] = _xcorr_max(corr_array, abs_max=False)

            misfit = corr_array[(len(corr_array) - 1) // 2 + int(shift_CC[iphase])]

            import numpy as np

            mid = (len(corr_array) - 1) / 2
            if len(corr_array) % 2 == 1:
                mid = int(mid)
            t = np.linspace(0, len(corr_array), len(corr_array), endpoint=False) - mid
            import matplotlib.pyplot as plt

            plt.close()
            plt.plot(t, corr_array)
            # plt.plot(st_syn[iphase].times(), st_syn[iphase].data)
            # plt.plot(st_obs[iphase].times(), st_obs[iphase].data)
            plt.show()
            plt.close()
            a = 1

        return misfit_CC


class Pol(_AbstractMisfit):
    name = "POL"
    description = "polarization based misfit"
    noise_level = False  # Do you need to determine level of noise for this misfit class?

    def __init__(self):
        pass

    def run_misfit(self, st_obs: _obspy.Stream, st_syn: _obspy.Stream):

        # TODO check this misfit
        P_Syn = np.zeros((3, len(st_synth.traces[0].data)))
        P_Syn[0, :] = st_synth.traces[0].data
        P_Syn[1, :] = st_synth.traces[3].data
        P_Syn[2, :] = st_synth.traces[5].data

        P_Obs = np.zeros((3, len(st_synth.traces[0].data)))
        P_Obs[0, :] = st_real.traces[0].data
        P_Obs[1, :] = st_real.traces[3].data
        P_Obs[2, :] = st_real.traces[5].data

        S_Syn = np.zeros((3, len(st_synth.traces[2].data)))
        S_Syn[0, :] = st_synth.traces[2].data
        S_Syn[1, :] = st_synth.traces[4].data
        S_Syn[2, :] = st_synth.traces[1].data

        S_Obs = np.zeros((3, len(st_synth.traces[2].data)))
        S_Obs[0, :] = st_real.traces[2].data
        S_Obs[1, :] = st_real.traces[4].data
        S_Obs[2, :] = st_real.traces[1].data

        Cov_P_syn = (1 / P_Syn.shape[1]) * np.matmul(P_Syn, P_Syn.T)
        Cov_P_obs = (1 / P_Obs.shape[1]) * np.matmul(P_Obs, P_Obs.T)
        Cov_S_syn = (1 / S_Syn.shape[1]) * np.matmul(S_Syn, S_Syn.T)
        Cov_S_obs = (1 / S_Obs.shape[1]) * np.matmul(S_Obs, S_Obs.T)

        from numpy import linalg as LA

        Eigval_Psyn, Eigvec_Psyn = LA.eig(Cov_P_syn)
        Eigval_Pobs, Eigvec_Pobs = LA.eig(Cov_P_obs)
        Eigval_Ssyn, Eigvec_Ssyn = LA.eig(Cov_S_syn)
        Eigval_Sobs, Eigvec_Sobs = LA.eig(Cov_S_obs)

        # TODO: check if you needed the maximum or the minimum eigenvalue
        max_ind_Eigval_Pobs = Eigval_Pobs.argmax()
        max_ind_Eigval_Psyn = Eigval_Pobs.argmax()

        # TODO: Use this value now as a misfit...
        delta_s = np.rad2deg(
            np.arccos(
                (
                    Eigvec_Pobs[:, max_ind_Eigval_Pobs][:, None].T
                    @ Eigvec_Psyn[:, max_ind_Eigval_Psyn][:, None]
                ).item()
            )
        )

