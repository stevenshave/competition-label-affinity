"""
Produce plot of competition experiment sensitivity, inhibitor KD x-axis

Generates a plot supporting the manuscript "Identification of optimum ligand affinity
for competition-based primary screens" by Shave et.al.
"""

from inspect import CO_ASYNC_GENERATOR
from claffinity import CompetitionLabelAffinity
cla=CompetitionLabelAffinity()
cla.plot_inhibitor_KD_vs_FLB()
