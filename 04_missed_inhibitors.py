"""
Calculate readout for high and low ligand affinity competition systems

Output the concentration and fraction ligand bound for two competition systems,
one using a high affinity (1 nM) ligand and the other low (100 nM).
Used in the preparation of the manuscript "Identification of optimum ligand
affinity for competition-based primary screens" by Shave et.al.
"""

import sys
from high_accuracy_binding_equations import *
from matplotlib import pyplot as plt
import numpy as np


# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

TARGET_FRACTION_BOUND = 0.7
LIGAND_CONC = 0.01
LOW_AFFINITY_LIGAND_KD = 0.100
HIGH_AFFINITY_LIGAND_KD = 0.001
INHIBITOR_KD=1.0
INHIBITOR_CONC = 10

low_affinity_ligand_system = {'kdpl': LOW_AFFINITY_LIGAND_KD, 'l': LIGAND_CONC, 'i': INHIBITOR_CONC, 'kdpi': INHIBITOR_KD}
low_affinity_ligand_system['p'] = calc_amount_p(TARGET_FRACTION_BOUND, LIGAND_CONC, low_affinity_ligand_system['kdpl'])
high_affinity_ligand_system = {'kdpl': HIGH_AFFINITY_LIGAND_KD,  'l': LIGAND_CONC, 'i': INHIBITOR_CONC, 'kdpi': INHIBITOR_KD}
high_affinity_ligand_system['p'] = calc_amount_p(
    TARGET_FRACTION_BOUND, LIGAND_CONC, high_affinity_ligand_system['kdpl'])

pl_low_affinity_ligand_system = competition_pl(**low_affinity_ligand_system)
pl_high_affinity_ligand_system = competition_pl(**high_affinity_ligand_system)

print(f"Low affinity ligand system: {low_affinity_ligand_system}, [PL]={pl_low_affinity_ligand_system:.4f}")
print(f"\tFracton ligand bound = {pl_low_affinity_ligand_system/low_affinity_ligand_system['l']}")
print()
print(f"High affinity ligand system: {high_affinity_ligand_system}, [PL]={pl_high_affinity_ligand_system:.4f}")
print(f"\tFracton ligand bound = {pl_high_affinity_ligand_system/high_affinity_ligand_system['l']}")
