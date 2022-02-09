"""
Produce plot of competition experiment sensitivity with fixed lignand

Generates a supporting plot of signal (PL) in a competition experiment vs
ligand KD for a fixed concentration and KD of inhibitor "Identification
of optimum ligand affinity for competition-based primary screens" by Shave et.al.
"""

from claffinity import CompetitionLabelAffinity
cla=CompetitionLabelAffinity()
cla.plot_signal_vs_ligand_kd_fixed_i()
