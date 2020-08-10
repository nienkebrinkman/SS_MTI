# from .GetData import GetData
# from .Inversion import Inversion
# from .PullData import PullData
# from .PhaseTracer import PhaseTracer
# from . import GreensFunctions
# from . import Forward
# from .PreProcess import PreProcess

# Dit is de manier om de module zelf te skippen en alleen de class te pakken
# from SS_MTI.PhaseTracer import PhaseTracer as PhaseTracer

# Dit pakt de hele module met onderliggende functions ZONDER _
from SS_MTI import (
    DataGetter,
    Forward,
    GreensFunctions,
    Inversion,
    PreProcess,
    PullData,
    Misfit,
    PostProcessing,
    MTDecompose,
    SourceTimeFunction,
    Read_H5,
)


__all__ = [
    "DataGetter",
    "Forward",
    "GreensFunctions",
    "Inversion",
    "PreProcess",
    "PullData",
    "PhaseTracer",
    "Misfit",
    "PostProcessing",
    "MTDecompose",
    "SourceTimeFunction",
    "Read_H5",
]
