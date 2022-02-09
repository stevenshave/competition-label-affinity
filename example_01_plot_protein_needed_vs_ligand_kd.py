"""
Produce plot of amount of protein needed over a range of ligand KDs

Generates a plot of protein needed to reach a desired fraction bound
in a competition experiment, dependant on ligand KD.  For the manuscript
"Identification of optimum ligand affinity for competition-based primary screens"
by Shave et.al.
"""

from claffinity import CompetitionLabelAffinity

competition_label_affinity=CompetitionLabelAffinity()
competition_label_affinity.plot_protein_needed_vs_ligand_kd()
