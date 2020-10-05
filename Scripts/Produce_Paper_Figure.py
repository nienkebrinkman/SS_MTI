"""
Script to generate a Figure for paper.
"""
__author__ = "Nienke Brinkman"

from os.path import join as pjoin
import matplotlib.pyplot as plt

# import matplotlib.image as mpimg
# from svgutils.compose import *

import SS_MTI
from SS_MTI import PostProcessing as _PostProcessing

save_folder = (
    "/home/nienke/Documents/Research/Data/MTI/Inversion/Result_2/5phases_cluster/Test_2020/"
)

depth = 38
fmin = 0.1
fmax = 0.4
event_name = "S0173a"
veloc_name = "TAYAK"
amount_of_phases = 5

""" Open misfit image """
misfit_file = pjoin(
    save_folder, f"Misfit_vs_Depth_{event_name}_{fmin}_{fmax}_L2_{veloc_name}.svg",
)

""" Open GS image """
GS_wf_file = pjoin(save_folder, f"GS_waveforms_{event_name}_{depth}_L2_{veloc_name}_Post.svg",)
GS_BBB_file = pjoin(save_folder, f"GS_BBB__{event_name}_{depth}_L2_{veloc_name}_Post.svg",)

""" Open Direct image """
Direct_wf_file = pjoin(
    save_folder, f"Direct_waveforms_{event_name}_{depth}_L2_{veloc_name}_Post.svg",
)
Direct_BBB_file = pjoin(save_folder, f"Direct_BB_{event_name}_{depth}_L2_{veloc_name}_Post.svg",)


import svgutils.transform as sg
import sys

# create new SVG figure
width = 65
height = 25
fig = sg.SVGFigure(f"{width}cm", f"{height}cm")

# load matpotlib-generated figures
move_X1 = (width * 0.32) / 0.026
move_X2 = (width * 0.3) / 0.026 + (width * 0.35) / 0.026


fig1 = sg.fromfile(misfit_file)
plot1 = fig1.getroot()
plot1.moveto(0, 4.1 / 0.026, scale=1.4)
txta = sg.TextElement(0, 1.3 / 0.026, "(a)", size=60, weight="bold")
txtb = sg.TextElement(move_X1, 1.3 / 0.026, "(b)", size=60, weight="bold")
txtc = sg.TextElement(move_X2, 1.3 / 0.026, "(c)", size=60, weight="bold")

txt1 = sg.TextElement(0 + 0.2 * move_X1, 2 / 0.026, f"Event {event_name}", size=60, weight="bold")

fig2 = sg.fromfile(GS_wf_file)
plot2 = fig2.getroot()
plot2.moveto(move_X1, 6.7 / 0.026, scale=1.2)

fig2a = sg.fromfile(GS_BBB_file)
plot2a = fig2a.getroot()
plot2a.moveto(move_X1 + 0.45 * move_X1, 1.0 / 0.026, scale=0.6)

txt2a = sg.TextElement(
    move_X1 + 0.45 * move_X1, 0.9 / 0.026, "Grid Search", size=40, weight="bold"
)
txt2b = sg.TextElement(
    move_X1 + 0.45 * move_X1, 2 / 0.026, f"Depth {depth}", size=35, weight="bold"
)


fig3 = sg.fromfile(Direct_wf_file)
plot3 = fig3.getroot()
plot3.moveto(move_X2, 6.7 / 0.026, scale=1.2)

fig3a = sg.fromfile(Direct_BBB_file)
plot3a = fig3a.getroot()
plot3a.moveto(move_X2 + 0.25 * move_X1, 1 / 0.026, scale=0.6)
txt3a = sg.TextElement(move_X2 + 0.6 * move_X1, 1 / 0.026, "Direct", size=40, weight="bold")
txt3b = sg.TextElement(
    move_X2 + 0.6 * move_X1, 2 / 0.026, f"Depth {depth}km", size=35, weight="bold"
)


fig.append([plot1, plot2, plot2a, plot3, plot3a])
fig.append([txt1, txt2a, txt2b, txt3a, txt3b, txta, txtb, txtc])
fig.save(pjoin(save_folder, f"{event_name}_{depth}_{veloc_name}.svg"))
# from IPython.display import SVG, display

# display(SVG(filename=pjoin(save_folder, "Test.svg")))

# here starts the assembling using svgutils
# Figure("16cm", "16cm", SVG(GS_wf_file)).save(pjoin(save_folder, "Test.svg"))
# SVG(pjoin(save_folder, "Test.svg"))

