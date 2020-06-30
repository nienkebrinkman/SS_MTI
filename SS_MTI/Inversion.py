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
from obspy.taup import TauPyModel
from obspy import UTCDateTime as utct


class Inversion:
    def __init__(
        self, forward_method: str, forward_dict:dict, rec_lat: float, rec_lon: float,
    ):
        """
        :param forward: string defining the forward modeller: "INSTASEIS", "REFLECTIVITY"
        :param forward_dict: Dict with all specification of the forward modeller (see input.toml)
        :param rec_lat: latitude of receiver station
        :param rec_lon: longitude of receiver station
        """
        if forward_method == "INSTASEIS":
            self.rec = instaseis.Receiver(latitude=rec_lat, longitude=rec_lon)
            self.veloc = instaseis.open_db(forward_dict["VELOC"])
            self.taup_veloc = TauPyModel(forward_dict["VELOC_taup"])
        elif forward_method == "REFLECTIVITY":
            raise ValueError("REFLECTIVITY CODE NEEDS TO BE IMPLEMENTED!!")
        else:
            raise ValueError(
                "forward_method can be either INSTASEIS or REFLECTIVITY in [FORWARD] of .toml file"
            )
        pass

    def Grid_Search(self, event: obspy.core.event.Event,depths: [float],strikes : [float], dips:[float], rakes: [float]):
        """
        Grid search over strike, dip, rake angles
        :param event: Obspy.event including waveforms and phase arrivals
        """

        for depth in depths:
            for strike in strikes:
                for dip in dips:
                    for rake in rakes:
                        pass



        pass

    def Direct(self, event: obspy.core.event.Event):
        pass

    def MH(self, event: obspy.core.event.Event):
        pass
