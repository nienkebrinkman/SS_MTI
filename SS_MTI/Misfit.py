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
        # sigmas_model = []
        for i in range(len(phases)):
            d_obs = _np.expand_dims(st_obs[i].data, axis=1)
            d_syn = _np.expand_dims(st_syn[i].data, axis=1)

            start_weight = self.weights[i][0]
            end_weight = self.weights[i][1]

            d_weight = _np.zeros_like(st_obs[i].data)
            d_weight[:samps] = start_weight
            d_weight[samps:] = end_weight

            # inv_Std = _np.diag(1 / (_np.std(d_obs) ** 2 * d_weight))
            # misfit_L2.append(0.5 * ((d_obs - d_syn).T @ inv_Std @ (d_obs - d_syn))[0][0])

            # TODO: determine sigma based on RMS of data + sigma noise!!
            # sigmas_model.append(_np.sqrt(_np.mean(d_syn ** 2 + d_obs ** 2)))

            misfit_L2.append(
                0.5
                * (
                    (d_obs - d_syn).T
                    @ (_np.expand_dims(1 / (sigmas[i] * d_weight), axis=1) * (d_obs - d_syn))
                )[0][0]
            )
            # misfit_L2.append(0.5 * ((d_obs - d_syn).T @ (d_obs - d_syn))[0][0])
        print(sigmas)
        print(misfit_L2)
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
    noise_level = True  # Do you need to determine level of noise for this misfit class?

    def __init__(
        self, components: [str], start_weight_len: float, weights: [_List[int]], dt: float
    ):
        self.components = components
        self.start_weight_len = start_weight_len
        self.weights = weights
        self.dt = dt

    def run_misfit(
        self, phases: [str], st_obs: _obspy.Stream, st_syn: _obspy.Stream, sigmas: [float]
    ):
        seen = set()
        uniq = [x for x in phases if x not in seen and not seen.add(x)]
        for i, uniq_phase in enumerate(uniq):
            inds_phase = [i for i, x in enumerate(phases) if x == uniq_phase]
            syn_unique = _obspy.Stream()
            obs_unique = _obspy.Stream()
            comp_order = []
            for ind, ind_phase in enumerate(inds_phase):
                comp_order.append(self.components[ind_phase])
                syn_unique += st_syn[ind_phase].copy()
                obs_unique += st_obs[ind_phase].copy()
            if i == 0:
                syn = _np.zeros((3, len(syn_unique.traces[0].data)))
                syn[0, :] = syn_unique.select(channel="BX" + "Z").copy()[0].data
                syn[1, :] = syn_unique.select(channel="BX" + "R").copy()[0].data
                syn[2, :] = syn_unique.select(channel="BX" + "T").copy()[0].data

                obs = _np.zeros((3, len(obs_unique.traces[0].data)))
                obs[0, :] = obs_unique.select(channel="BH" + "Z").copy()[0].data
                obs[1, :] = obs_unique.select(channel="BH" + "R").copy()[0].data
                obs[2, :] = obs_unique.select(channel="BH" + "T").copy()[0].data
            else:
                syn_temp = _np.zeros((3, len(syn_unique.traces[0].data)))
                syn_temp[0, :] = syn_unique.select(channel="BX" + "Z").copy()[0].data
                syn_temp[1, :] = syn_unique.select(channel="BX" + "R").copy()[0].data
                syn_temp[2, :] = syn_unique.select(channel="BX" + "T").copy()[0].data

                syn = _np.hstack((syn, syn_temp))

                obs_temp = _np.zeros((3, len(obs_unique.traces[0].data)))
                obs_temp[0, :] = obs_unique.select(channel="BH" + "Z").copy()[0].data
                obs_temp[1, :] = obs_unique.select(channel="BH" + "R").copy()[0].data
                obs_temp[2, :] = obs_unique.select(channel="BH" + "T").copy()[0].data

                obs = _np.hstack((obs, obs_temp))

        Cov_syn = (1 / syn.shape[1]) * (syn @ syn.T)
        Cov_obs = (1 / obs.shape[1]) * (obs @ obs.T)

        from numpy import linalg as LA

        Eigval_syn, Eigvec_syn = LA.eig(Cov_syn)
        Eigval_obs, Eigvec_obs = LA.eig(Cov_obs)

        # TODO: check if you needed the maximum or the minimum eigenvalue
        max_ind_Eigval_obs = Eigval_obs.argmax()
        max_ind_Eigval_syn = Eigval_syn.argmax()

        # TODO: Use this value now as a misfit...
        delta_s = _np.rad2deg(
            _np.arccos(
                (
                    Eigvec_obs[:, max_ind_Eigval_obs][:, None].T
                    @ Eigvec_syn[:, max_ind_Eigval_syn][:, None]
                ).item()
            )
        )
        L2 = [0] * len(phases)
        L2[0] = delta_s
        return L2
