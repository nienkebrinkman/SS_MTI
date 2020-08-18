""" This class is based on TauP calculator"""

from obspy.taup import TauPyModel as _TauPyModel


def get_traveltime(
    model: _TauPyModel, phase: str, depth: float, distance: float, take_off: bool = False
) -> float:
    """
    Get travel time of phase
    :param model: Velocity model
    :param phases: list of phases to include in the inversion
    :param depth: Depth of event in km
    :param distance: Distance of event in degrees
    :param take_off: return take-off angle instead of traveltime
    """
    try:
        tt = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=distance, phase_list=[phase]
        )
        if take_off:
            return tt[0].takeoff_angle
        else:
            return tt[0].time  # TODO: CHANGE THIS
    except IndexError as e:
        # raise e("{} not arriving at {}km depth and {} degrees".format(phase, depth, distance))
        print("{} not arriving at {}km depth and {} degrees".format(phase, depth, distance))
        return None

