import abc
import instaseis
from obspy.taup import TauPyModel as _TauPyModel
import obspy
from typing import Union as _Union, Tuple as _Tuple

from SS_MTI import PhaseTracer as _PhaseTracer
from SS_MTI import GreensFunctions as _GreensFunctions


class _AbstractForward(metaclass=abc.ABCMeta):

    name = "abstract forward model"

    def __init__(self):
        pass

    @abc.abstractmethod
    def greens_functions(self):
        pass

    @abc.abstractmethod
    def generate_synthetic_data(self):
        pass


class Instaseis(_AbstractForward):

    name = "Instaseis based forward model"

    def __init__(
        self,
        instaseis_db: instaseis.open_db,
        taup_model: _TauPyModel,
        rec_lat: float,
        rec_lon: float,
        or_time: obspy.UTCDateTime,
        dt: float = 0.05,
    ) -> None:
        """ Setup of instaseis forward modeling is specified """
        self.rec = instaseis.Receiver(latitude=rec_lat, longitude=rec_lon)
        self.veloc = instaseis.open_db(instaseis_db)
        self.taup_veloc = _TauPyModel(taup_model)
        self.or_time = or_time
        self.dt = dt

    def greens_functions(
        self,
        phase: str,
        comp: str,
        depth: float,
        distance: float,
        lat_src: float,
        lon_src: float,
        tstar: _Union[float, str],
        LQT: bool = False,
        inc: float = None,
        baz: float = None,
    ) -> _Tuple[obspy.Stream, float]:
        """ GREENS FUNCTIONS FOR INSTASEIS FORWARD MODELING: """

        """ SYNTHETIC TRAVEL TIME: """
        syn_tt = _PhaseTracer.get_traveltime(self.taup_veloc, phase, depth, distance)

        """ GREEN'S FUNCTIONS: """
        st_GF = _GreensFunctions.make_GF(
            or_time=self.or_time,
            lat_src=lat_src,
            lon_src=lon_src,
            depth=depth,
            distance=distance,
            rec=self.rec,
            db=self.db,
            dt=self.dt,
            comp=comp,
            tstar=tstar,
            LQT=LQT,
            inc=inc,
            baz=baz,
        )
        return st_GF, syn_tt

    def generate_synthetic_data(self):
        st_syn = obspy.Stream

        for i, phase in enumerate(phases):

            syn_tr = self.GF.from_GF(syn_GF[i], strike, dip, rake)
            phase_tr = syn.slice(
                starttime=event.origin_time + syn_tt[i] - t_pre[iphase],
                endtime=event.origin_time + syn_tt[i] + t_post[iphase],
            )


class reflectivity(_AbstractForward):

    name = "reflectivity based forward model"

    def __init__(self):
        raise ValueError("REFLECTIVITY CODE NEEDS TO BE IMPLEMENTED!!")

    def greens_functions(self, depth: float, phases: [str]):
        """ GREENS FUNCTIONS FOR REFLECTIVITY FORWARD MODELING: """
        pass

    def generate_synthetic_data(self):
        """ GENERATING SYNTHETIC WAVEFORMS: """
        pass
