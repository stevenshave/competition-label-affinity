"""
Produce plot of competition experiment sensitivity, ligand KD x-axis

Generates a plot used in the manuscript "Revisiting labelled ligand affinity
in competition experiments" by Shave et.al.
"""

import sys
from matplotlib import pyplot as plt
import numpy as np
from high_accuracy_binding_equations import *


# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 2000 # Publication used 2000 pts along X
TARGET_FRACTION_L_BOUND = 0.7

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10**(-x_axis)  
inhibitor_kds = np.array([1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6])
protein_concs = np.array(calc_amount_p(
    TARGET_FRACTION_L_BOUND, 10e-9, ligand_kd_range))
y = np.full((inhibitor_kds.shape[0], x_axis.shape[0]), np.nan)

for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
    for it_ligand_kd_range, ligand_kd in enumerate(ligand_kd_range):
        lig_conc = 10e-9
        i_conc = 10e-6
        p = protein_concs[it_ligand_kd_range]
        y[it_inhibitor_kds][it_ligand_kd_range] = competition_pl(
            **{'p': p, 'l': lig_conc, 'i': i_conc, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/lig_conc


fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

plot_line_labels = [
    r'K$_\mathrm{D}$PI=1 nM',
    r'K$_\mathrm{D}$PI=10 nM',
    r'K$_\mathrm{D}$PI=100 nM',
    r'K$_\mathrm{D}$PI=1 µM',
    r'K$_\mathrm{D}$PI=10 µM',
    r'K$_\mathrm{D}$PI=100 µM',
]

plot_marker_styles = list(reversed([
    '',
    '+',
    'x',
    'd',
    "s",
    "8",
]))

for i in reversed(range(inhibitor_kds.shape[0])):
    ax.plot(x_axis, y[i], 'k', label=plot_line_labels[i],
            marker=plot_marker_styles[i], markevery=NUM_POINTS_ON_XAXIS//20, linewidth=1)
ax.set_xlabel(r"Ligand pK$_\mathrm{D}$")
ax.set_ylabel("Fraction ligand bound")
ax.legend()
ax.grid()
ax.title.set_text(r"Protein-ligand signal over a range of ligand K$_\mathrm{D}$s, [L$_0$]=10 nM, [I$_0$]=10 " +
                  r"$\mathrm{\mu}$M"+f"\nTarget fraction ligand bound without inhibitor = {TARGET_FRACTION_L_BOUND}")
ax.set_xlim(3, 12)
ax.set_ylim(0, TARGET_FRACTION_L_BOUND*1.1)
plt.show()
