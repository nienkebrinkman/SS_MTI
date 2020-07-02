import obspy


def filter(tr, fmin=1.0 / 10.0, fmax=1.0 / 2, zerophase=True):
    tr.filter("highpass", freq=fmin, zerophase=zerophase)
    tr.filter("highpass", freq=fmin, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, zerophase=zerophase)
    tr.filter("lowpass", freq=fmax, zerophase=zerophase)


def prepare_event_data(
    event: obspy.core.event.Event,
    phases: [str],
    components: [str],
    tts: [float],
    t_pre: [float],
    t_post: [float],
    fmin: float = 1.0 / 8,
    fmax: float = 1.0 / 3.5,
    zerophase: bool = False,
) -> obspy.Stream:
    """
    Pre-process data that is specified in the event object.
    """
    st_obs = obspy.Stream()
    for i, phase in enumerate(phases):
        if components[i] in ["Z", "N", "E"]:
            tr_orig = event.waveforms_VBB.select(channel="BH" + components[i]).copy()[0]
        else:
            st_raw = event.waveforms_VBB.copy()
            st_raw.rotate(method="NE->RT", back_azimuth=event.baz)
            tr_orig = st_raw.select(channel="BH" + components[i])[0]

        filter(tr_orig, fmin=fmin, fmax=fmax, zerophase=zerophase)

        tr_window = tr_orig.slice(
            starttime=event.origin_time + tts[i] - t_pre[i],
            endtime=event.origin_time + tts[i] + t_post[i],
        )

        st_obs += tr_window
    return st_obs
