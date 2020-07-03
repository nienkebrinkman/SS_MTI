import abc
import instaseis
from obspy.taup import TauPyModel as _TauPyModel
import obspy
from typing import Union as _Union, Tuple as _Tuple

from SS_MTI import PhaseTracer as _PhaseTracer
from SS_MTI import GreensFunctions as _GreensFunctions
from SS_MTI import PreProcess as _PreProcess


class _AbstractForward(metaclass=abc.ABCMeta):

    name = "abstract forward model"

    def __init__(self):
        pass

    @abc.abstractmethod
    def get_greens_functions(self):
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
        start_cut: float = 100.0,
        end_cut: float = 800.0,
    ) -> None:
        """ Setup of instaseis forward modeling is specified """
        self.rec = instaseis.Receiver(latitude=rec_lat, longitude=rec_lon)
        self.db = instaseis.open_db(instaseis_db)
        self.taup_veloc = _TauPyModel(taup_model)
        self.or_time = or_time
        self.dt = dt
        self.start_cut = start_cut
        self.end_cut = end_cut

    def get_greens_functions(
        self,
        phase: str,
        comp: str,
        depth: float,
        distance: float,
        lat_src: float,
        lon_src: float,
        tstar: _Union[float, str],
        M0: float,
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
            M0=M0,
        )
        return st_GF, syn_tt

    def generate_synthetic_data(
        self,
        st_GF: obspy.Stream,
        focal_mech: [float],
        M0: float,
        slice: bool = False,
        tt: float = None,
        t_pre: float = None,
        t_post: float = None,
        filter: bool = False,
        fmin: float = None,
        fmax: float = None,
        zerophase: bool = False,
    ):
        """ Generate synthetic waveforms 
        :param st_GF: 
        :param focal_mech: strike,dip,rake or m_rr, m_pp, m_tt, m_rp, m_rt, m_tp
        :param M0: scalar moment
        :param slice: if true the trace will be slices around the tt
        :param tt: travel time in seconds after origin time
        """

        assert (
            slice == True and tt is not None and t_pre is not None and t_post is not None
        ), "if slice is set to True you have to specify tt, t_pre and t_post"

        syn_tr_full = _GreensFunctions.from_GF(st_GF, focal_mech, M0)
        syn_tr = syn_tr_full.slice(
            starttime=self.or_time + tt - t_pre, endtime=self.or_time + tt + t_post,
        )

        if filter:
            assert (
                fmin is not None and fmax is not None
            ), "if filter == True, specify fmin, fmax and zerophase"
            _PreProcess.filter_tr(syn_tr, fmin=fmin, fmax=fmax, zerophase=zerophase)

        return syn_tr


class reflectivity(_AbstractForward):

    name = "reflectivity based forward model"

    def __init__(self):
        raise ValueError("REFLECTIVITY CODE NEEDS TO BE IMPLEMENTED!!")

    def get_greens_functions(self, depth: float, phases: [str]):
        """ GREENS FUNCTIONS FOR REFLECTIVITY FORWARD MODELING: """
        pass

    def generate_synthetic_data(self):
        """ GENERATING SYNTHETIC WAVEFORMS: """
        pass
