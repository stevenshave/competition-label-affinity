"""
Produce plot of competition experiment sensitivity, ligand KD x-axis

Generates a plot used in the manuscript "Identification of optimum ligand affinity
for competition-based primary screens" by Shave et.al.
"""
from claffinity import CompetitionLabelAffinity
cla=CompetitionLabelAffinity()
cla.plot_ligand_KD_vs_FLB()
