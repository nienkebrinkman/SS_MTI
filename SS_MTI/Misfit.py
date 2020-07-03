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
    name = "L2 based misfit"
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
        self, phases: [str], st_obs: _obspy.Stream, st_syn: _obspy.Stream, sigmas: [float]
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
        for i in range(len(st_obs)):
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
                )
            )
        return misfit_L2


class CC(_AbstractMisfit):
    name = "cross-correlation based misfit"
    noise_level = False  # Do you need to determine level of noise for this misfit class?

    def __init__(self):
        pass

    def run_misfit(self):
        pass


class Pol(_AbstractMisfit):
    name = "polarization based misfit"
    noise_level = False  # Do you need to determine level of noise for this misfit class?

    def __init__(self):
        pass

    def run_misfit(self):
        pass

