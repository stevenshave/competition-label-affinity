"""
Produce plot of competition experiment sensitivity with fixed lignand

Generates a supporting plot of signal (PL) in a competition experiment vs
ligand KD for a fixed concentration and KD of inhibitor "Revisiting labelled
ligand affinity in competition experiments" by Shave et.al.
"""

import sys
from matplotlib import pyplot as plt
import numpy as np
from high_accuracy_binding_equations import *


# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 2000 # Publication used 2000 pts along X
LIGAND_CONC=10e-9
INHIBITOR_CONC=10e-6

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10**(-x_axis)
protein_conc = 1e-6
y = np.full((x_axis.shape[0]), np.nan)
for i in range(ligand_kd_range.shape[0]):
    y[i] = competition_pl(**{'p': protein_conc, 'l': LIGAND_CONC, 'i': INHIBITOR_CONC, 'kdpl': ligand_kd_range[i], 'kdpi': 1e-6})/LIGAND_CONC

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

ax.plot(x_axis, y, 'k')
ax.set_xlabel(r"Ligand pK$_\mathrm{D}$")
ax.set_ylabel("Fraction ligand bound")
ax.title.set_text(r"Competition experiment fraction ligand bound over a range of ligand K$_\mathrm{D}$s"+"\n"+r"[P]=10 µM, [L]=10 nM, [I]=10 µM, inhibitor K$_\mathrm{D}$=10 µM")
ax.set_xlim(3, 12)
ax.grid()
ax.set_ylim(0, 1.001)
plt.show()
