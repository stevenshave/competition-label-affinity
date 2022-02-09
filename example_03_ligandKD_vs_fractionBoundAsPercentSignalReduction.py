"""
Produce plot of competition experiment sensitivity in terms of % signal

Generates a plot which was explored but not used in the preparation of the
manuscript "Identification of optimum ligand affinity for competition-based
primary screens" by Shave et.al.
Expresses PL complex in terms of maximum achievable signal.
"""


from claffinity import CompetitionLabelAffinity
cla=CompetitionLabelAffinity()
cla.plot_ligand_kd_vs_FLB_as_percentage()
