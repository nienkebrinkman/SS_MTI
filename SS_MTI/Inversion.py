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
from typing import List as _List, Union as _Union


from SS_MTI import PreProcess as _PreProcess
from SS_MTI import Forward as _Forward
from SS_MTI import Misfit as _Misfit


# def Grid_Search_run(fwd: _Forward._AbstractForward):

#     print(f"Running grid search with model: {fwd.name}")

#     _preprocess()

#     "Actual algorithm"

#     for i in grid:
#         fwd.greens_functions()

#     _postprocess()


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
):
    """
    Grid search over strike, dip, rake angles
    :param event: Obspy.event including waveforms and phase arrivals
    :param phases: list of phases to include in the inversion
    """
    print(f"Running grid search with model: {fwd.name}")
    print(f"and with misfit: {misfit.name}")
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

    # TODO: IMPLEMENT LQT COORDINATE SYSTEM
    LQT_value = False
    baz = None
    inc = None

    """ step 1: PRE-PROCESS THE OBSERVED DATA """
    obs_tt = []
    for i, phase in enumerate(phases):
        obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])

    st_obs, sigmas = _PreProcess.prepare_event_data(
        event=event,
        phases=phases,
        components=components,
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
        """ Step 2: GENERATE GREEN'S FUNCTION AT SPECIFIC DEPTH """
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
            )
            # TODO: filter the streams
            syn_GFs.append(syn_GF)
            syn_tts.append(syn_tt)

        for strike in strikes:
            for dip in dips:
                for rake in rakes:

                    focal_mech = [strike, dip, rake]
                    st_syn = obspy.Stream()

                    for i, phase in enumerate(phases):
                        tr_syn = fwd.generate_synthetic_data(
                            st_GF=syn_GFs[i],
                            focal_mech=focal_mech,
                            M0=M0,
                            slice=True,
                            tt=syn_tts[i],
                            t_pre=t_pre[i],
                            t_post=t_post[i],
                            filter=filter_par,
                            fmin=fmin,
                            fmax=fmax,
                        )
                        st_syn += tr_syn

                    chi = misfit.run_misfit(
                        phases=phases, st_obs=st_obs, st_syn=st_syn, sigmas=sigmas
                    )

                    a = 1


def Direct(self, event: obspy.core.event.Event):
    pass


def MH(self, event: obspy.core.event.Event):
    pass


class Inversion:
    def __init__(
        self, forward_method: str, forward_dict: dict, rec_lat: float, rec_lon: float,
    ):
        """
        :param forward: string defining the forward modeller: "INSTASEIS", "REFLECTIVITY"
        :param forward_dict: Dict with all specification of the forward modeller (see input.toml)
        :param rec_lat: latitude of receiver station
        :param rec_lon: longitude of receiver station
        """
        if forward_method == "INSTASEIS":
            self.forward = _Forward.Instaseis(
                instaseis_db=forward_dict["VELOC"],
                taup_model=forward_dict["VELOC_taup"],
                rec_lat=rec_lat,
                rec_lon=rec_lon,
            )
        elif forward_method == "REFLECTIVITY":
            self.forward = _Forward.reflectivity()
        else:
            raise ValueError(
                "forward_method can be either INSTASEIS or REFLECTIVITY in [FORWARD] of .toml file"
            )
        pass

    def Grid_Search(
        self,
        event: obspy.core.event.Event,
        phases: [str],
        phase_corrs: [float],
        components: [str],
        t_pre: [str],
        t_post: [str],
        depths: [float],
        strikes: [float],
        dips: [float],
        rakes: [float],
        tstars: _Union[_List[float], _List[str]] = None,
        fmin: float = 1.0 / 8,
        fmax: float = 1.0 / 3.5,
        zerophase: bool = False,
    ):
        """
        Grid search over strike, dip, rake angles
        :param event: Obspy.event including waveforms and phase arrivals
        :param phases: list of phases to include in the inversion
        """

        if tstars is None:
            tstars = [None] * len(phases)

        # TODO: IMPLEMENT LQT COORDINATE SYSTEM
        LQT_value = False

        """ step 1: PRE-PROCESS THE OBSERVED DATA """
        proc = _PreProcess()

        obs_tt = []
        for i, phase in enumerate(phases):
            obs_tt.append(utct(event.picks[phase]) - event.origin_time + phase_corrs[i])
        st_obs = proc.prepare_event_data(
            event=event,
            phases=phases,
            components=components,
            tts=obs_tt,
            t_pre=t_pre,
            t_post=t_post,
            fmin=fmin,
            fmax=fmax,
        )

        for depth in depths:
            """ Step 2: calculate synthetic arrivals """

            self.forward.greens_functions()

            for strike in strikes:
                for dip in dips:
                    for rake in rakes:

                        pass

        pass

    def Direct(self, event: obspy.core.event.Event):
        pass

    def MH(self, event: obspy.core.event.Event):
        pass

