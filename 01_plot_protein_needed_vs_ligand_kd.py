"""
Produce plot of amount of protein needed over a range of ligand KDs

Generates a plot of protein needed to reach a desired fraction bound
in a competition experiment, dependant on ligand KD.  For the manuscript
"Revisiting labelled ligand affinity in competition experiments" by Shave 
et.al.
"""

import sys
from matplotlib import pyplot as plt
import numpy as np
from high_accuracy_binding_equations import *

# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 1000 # Publication used 2000 pts along X
LIGAND_CONC=10e-9
TARGET_FRACTION_L_BOUND = 0.7

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10**(-x_axis)
print(ligand_kd_range)
protein_concs = np.array(calc_amount_p(TARGET_FRACTION_L_BOUND, LIGAND_CONC, ligand_kd_range))

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))
ax.plot(x_axis, protein_concs, 'k')
ax.set_xlabel(r"Ligand pK$_\mathrm{D}$")
ax.set_ylabel(r"[P$_0$] (M)")
plt.yscale("log")
ax.grid()
ax.title.set_text(r"Protein required for 0.7 fraction ligand bound vs labeled ligand K$_\mathrm{D}$, [L$_0$]=10 nM")
ax.set_xlim(3, 12)
plt.show()
