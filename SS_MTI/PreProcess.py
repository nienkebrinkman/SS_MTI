import obspy
import numpy as _np
from matplotlib import mlab as _mlab
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
from typing import Tuple as _Tuple


def Get_location(la_s, lo_s, la_r, lo_r, radius=3389.5, flattening=0):
    dist, az, baz = gps2dist_azimuth(
        lat1=la_s, lon1=lo_s, lat2=la_r, lon2=lo_r, a=radius, f=flattening
    )
    epi = kilometer2degrees(dist, radius=radius)
    return epi, az, baz


def calc_PSD(tr, winlen_sec):
    # tr.taper(0.05)
    Fs = tr.stats.sampling_rate
    winlen = min(winlen_sec * Fs, (tr.stats.endtime - tr.stats.starttime) * Fs / 2.0)
    NFFT = obspy.signal.util.next_pow_2(winlen)
    pad_to = _np.max((NFFT * 2, 1024))
    p, f = _mlab.psd(
        tr.data, Fs=Fs, NFFT=NFFT, detrend="linear", pad_to=pad_to, noverlap=NFFT // 2
    )
    return f, p


def filter_tr(tr, fmin=1.0 / 10.0, fmax=1.0 / 2, zerophase=False):
    tr.filter("highpass", freq=fmin, corners=4, zerophase=zerophase)
    tr.filter("highpass", freq=fmin, corners=4, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, corners=4, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, corners=4, zerophase=zerophase)


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
                starttime=event.origin_time + tts[i] - 15.0,
                endtime=event.origin_time + tts[i] - 5.0,
            )
            sigmas.append(_np.std(tr_noise.data))

            # Path = "/home/nienke/Documents/Research/Data/Noise/"
            # File_names = [
            #     "XB.02.ELYSE.BHE-2019.274T0809-2019.274T0920",
            #     "XB.02.ELYSE.BHN-2019.274T0809-2019.274T0920",
            #     "XB.02.ELYSE.BHZ-2019.274T0809-2019.274T0920",
            # ]
            # st_noise = obspy.Stream()

            # for file in File_names:
            #     tr = obspy.read(Path + file)
            #     st_noise += tr

            # if components == "LQT":
            #     raise ValueError("LQT orientation Not implemented yet")
            #     # TODO: implement LQT orientation
            # else:
            #     st_noise.rotate(method="NE->RT", back_azimuth=event.baz)

            # for trace in st_obs:

            #     chan = trace.stats.channel
            #     desired_dt = trace.stats.delta
            #     desired_npts = len(trace.data)
            #     noise_trace = st_noise.select(channel=chan)[0]
            #     noise = noise_trace.data[
            #         int(1200 / desired_dt) : int(1200 / desired_dt) + desired_npts
            #     ]
            #     sigmas.append(_np.std(noise))

        else:
            sigmas.append(None)

        st_obs += tr_window
    return st_obs, sigmas
