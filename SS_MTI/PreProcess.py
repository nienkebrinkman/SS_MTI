import obspy
import numpy as _np
from typing import Tuple as _Tuple


def filter_tr(tr, fmin=1.0 / 10.0, fmax=1.0 / 2, zerophase=False):
    tr.filter("highpass", freq=fmin, zerophase=zerophase)
    tr.filter("highpass", freq=fmin, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, zerophase=zerophase)


def rotate_event_data(
    event: obspy.core.event.Event, component: str,
):
    """ 
    Selects raw data from the event and rotates into component

    """
    if component in ["Z", "N", "E"]:
        tr_orig = event.waveforms_VBB.select(channel="BH" + component).copy()[0]
    else:
        st_raw = event.waveforms_VBB.copy()
        st_raw.rotate(method="NE->RT", back_azimuth=event.baz)
        tr_orig = st_raw.select(channel="BH" + component)[0]
    return tr_orig


def prepare_event_data(
    event: obspy.core.event.Event,
    phases: [str],
    components: [str],
    slice: bool = False,
    tts: [float] = None,
    t_pre: [float] = None,
    t_post: [float] = None,
    filter: bool = False,
    fmin: float = None,
    fmax: float = None,
    zerophase: bool = False,
    noise_level=False,
) -> _Tuple[obspy.Stream, obspy.Stream]:
    """
    Raw data is selected from event and rotated into components.
    The data can be filtered if specified.
    Then, the raw data is sliced around arrival times (tts) with t_pre and t_post
    Optional: slice the data before the arrival to obtain information about noise level,
    then you will have another output

    """
    if filter:
        assert (
            fmin is not None and fmax is not None
        ), "if filter == True, specify fmin, fmax and zerophase"

    assert (
        slice == True and tts is not None and t_pre is not None and t_post is not None
    ) or slice == False, "if slice is set to True you have to specify tt, t_pre and t_post"

    st_obs = obspy.Stream()
    if noise_level:
        print("noise stream will be a 30 second" " window starting 10 seconds before arrival")
    sigmas = []
    for i, phase in enumerate(phases):
        tr_orig = rotate_event_data(event, components[i])
        if filter:
            filter_tr(tr_orig, fmin=fmin, fmax=fmax, zerophase=zerophase)

        if slice:
            tr_window = tr_orig.slice(
                starttime=event.origin_time + tts[i] - t_pre[i],
                endtime=event.origin_time + tts[i] + t_post[i],
            )
        else:
            tr_window = tr_orig.slice(starttime=event.origin_time, endtime=tr_orig.stats.endtime,)
            # tr_window = tr_orig.copy()

        if noise_level:
            tr_noise = tr_orig.slice(
                starttime=event.origin_time + tts[i] - t_pre[i] - 40.0,
                endtime=event.origin_time + tts[i] - t_pre[i] - 10.0,
            )
            sigmas.append(_np.std(tr_noise.data))
        else:
            sigmas.append(None)

        st_obs += tr_window
    return st_obs, sigmas
