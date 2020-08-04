from obspy.core.event import Event
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from obspy import UTCDateTime as utct
from obspy.taup import TauPyModel
import obspy
import instaseis

from SS_MTI import PhaseTracer
from SS_MTI import GreensFunctions


class EventObj:
    @staticmethod
    def Get_location(la_s, lo_s, la_r, lo_r, radius=3389.5, flattening=0):
        dist, az, baz = gps2dist_azimuth(
            lat1=la_s, lon1=lo_s, lat2=la_r, lon2=lo_r, a=radius, f=flattening
        )
        epi = kilometer2degrees(dist, radius=radius)
        return epi, az, baz

    def __init__(
        self,
        or_time: utct = utct("2020-3-10T12:00:00"),
        lat_src: float = 10.99032013,
        lon_src: float = 170,
        lat_rec: float = 4.502384,
        lon_rec: float = 135.623447,
        depth: float = 45.0,
        name: str = "Test_Event",
    ):
        """
    Create a seismic event
    :param rec_lat: latitude receiver
    :param rec_lon: longitude receiver
    """
        self.event = Event()
        self.event.latitude = lat_src
        self.event.longitude = lon_src
        self.event.depth = depth
        self.event.name = name

        self.lat_rec = lat_rec
        self.lon_rec = lon_rec

        epi, az, baz = EventObj.Get_location(
            self.event.latitude, self.event.longitude, self.lat_rec, self.lon_rec
        )

        self.event.distance = epi
        print(self.event.distance)
        self.event.az = az
        self.event.baz = baz
        self.event.origin_time = or_time

    def add_picks(self, taup_model: TauPyModel, depth: float, phases: [str] = ["P", "S"]):
        self.event.picks = {}
        for phase in phases:
            self.event.picks[phase] = utct(
                self.event.origin_time
                + PhaseTracer.get_traveltime(
                    model=taup_model, phase=phase, depth=depth, distance=self.event.distance
                )
            )

    def add_waveforms(
        self,
        instaseis_db_path: str,
        focal_mech: [float],
        M0: float = None,
        dt: float = 0.05,
        components: str = "ZRT",
        kind: str = "displacement",
        noise: bool = False,
    ):
        """
    Add waveforms to event object using instaseis
    :param instaseis_db_path: path or url to instaseis database
    :param rec_lat: latitude receiver
    :param rec_lon: longitude receiver
    :param focal_mech: strike,dip,rake or m_rr, m_pp, m_tt, m_rp, m_rt, m_tp
    :param M0: scalar moment, only necessesary when focal_mech strike,dip,rake
    :param components: components of the seismogram (ZRT, ZNE, LQT)
    :param noise: real Martian noise will be added to the seismogram  
    """

        assert (M0 is None and len(focal_mech) == 6) or (
            M0 is not None and len(focal_mech) == 3
        ), (
            "focal_mech length is incorrect. "
            "If you specify M0, focal_mech is [strike,dip,rake]. "
            "Otherwise focal_mech is [m_rr, m_pp, m_tt, m_rp, m_rt, m_tp]"
        )

        receiver = instaseis.Receiver(
            latitude=self.lat_rec,
            longitude=self.lon_rec,
            network="XB",
            station="ELYSE",
            location="02",
        )

        db = instaseis.open_db(instaseis_db_path)

        if len(focal_mech) == 3:
            focal_mech = GreensFunctions.convert_SDR(
                focal_mech[0], focal_mech[1], focal_mech[2], M0
            )

        m_rr = focal_mech[0]
        m_pp = focal_mech[1]
        m_tt = focal_mech[2]
        m_rp = focal_mech[3]
        m_rt = focal_mech[4]
        m_tp = focal_mech[5]
        print(m_rr, m_pp, m_tt, m_rp, m_rt, m_tp)

        src = instaseis.Source(
            latitude=self.event.latitude,
            longitude=self.event.longitude,
            depth_in_m=self.event.depth * 1000,
            m_rr=m_rr,
            m_tt=m_tt,
            m_pp=m_pp,
            m_rt=m_rt,
            m_rp=m_rp,
            m_tp=m_tp,
            time_shift=None,
            sliprate=None,
            dt=None,
            origin_time=self.event.origin_time,
        )

        if components == "LQT":
            st_obs = db.get_seismograms(
                source=src, receiver=receiver, components=components, kind=kind, dt=dt
            )
            st_obs.rotate(method="RT->NE", back_azimuth=self.event.baz)
        else:
            st_obs = db.get_seismograms(
                source=src, receiver=receiver, components=components, kind=kind, dt=dt
            )

        st_obs[0].stats.channel = st_obs[0].stats.channel.replace("X", "H")
        st_obs[1].stats.channel = st_obs[1].stats.channel.replace("X", "H")
        st_obs[2].stats.channel = st_obs[2].stats.channel.replace("X", "H")

        st_obs.trim(starttime=self.event.origin_time, endtime=self.event.origin_time + 800.0)

        if noise:
            Path = "/home/nienke/Documents/Research/Data/Noise/"
            File_names = [
                "XB.02.ELYSE.BHE-2019.274T0809-2019.274T0920",
                "XB.02.ELYSE.BHN-2019.274T0809-2019.274T0920",
                "XB.02.ELYSE.BHZ-2019.274T0809-2019.274T0920",
            ]
            st_noise = obspy.Stream()

            for file in File_names:
                tr = obspy.read(Path + file)
                st_noise += tr

            if components == "LQT":
                raise ValueError("LQT orientation Not implemented yet")
                # TODO: implement LQT orientation
            else:
                st_noise.rotate(method="NE->RT", back_azimuth=self.event.baz)

            for trace in st_obs:
                chan = trace.stats.channel
                desired_dt = trace.stats.delta
                desired_npts = len(trace.data)
                noise_trace = st_noise.select(channel=chan)[0]
                noise = noise_trace.data[int(1200 / dt) : int(1200 / dt) + desired_npts]
                trace.data += noise

        self.event.waveforms_VBB = st_obs

