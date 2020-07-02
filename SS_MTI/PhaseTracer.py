""" This class is based on TauP calculator"""

from obspy.taup import TauPyModel as _TauPyModel


def get_traveltime(model: _TauPyModel, phase: str, depth: float, distance: float) -> float:
    """
    Get travel time of phase
    :param model: Velocity model
    :param phases: list of phases to include in the inversion
    :param depth: Depth of event in km
    :param distance: Distance of event in degrees
    """
    try:
        tt = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=distance, phase_list=[phase]
        )
    except IndexError as e:
        raise e(
            "{} not arriving at {}km depth and {} degrees".format(phase, depth, event.distance)
        )
    return tt[0].time

