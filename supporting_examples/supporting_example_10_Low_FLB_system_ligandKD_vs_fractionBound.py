"""
Produce plot of competition experiment sensitivity, inhibitor KD x-axis

Generates a plot used as supporting information Figure S8 in the manuscript
"Identification of optimum ligand affinity for competition-based primary screens"
by Shave et.al.
"""

import sys
from matplotlib import pyplot as plt
import numpy as np
from claffinity.high_accuracy_binding_equations import *


# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 2000 # Publication used 2000 pts along X
TARGET_FRACTION_L_BOUND = 0.7
LIGAND_CONC=100e-9
PROTEIN_CONC=8e-9

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
inhibitor_kd_range = 10**(-x_axis)  
ligand_kds = np.array([1e-6])
y = np.full((ligand_kds.shape[0], x_axis.shape[0]), np.nan)

for it_ligand_kd, ligand_kd in enumerate(ligand_kds):
    print(f"Generating: {it_ligand_kd+1}/{len(ligand_kds)}")
    for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kd_range):
        i_conc = 10e-6
        y[it_ligand_kd][it_inhibitor_kds] = competition_pl(
            **{'p': PROTEIN_CONC, 'l': LIGAND_CONC, 'i': i_conc, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/LIGAND_CONC

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

plot_line_labels = [
    r'K$_\mathrm{D}$PL=1 µM',
]

plot_marker_styles = list(reversed([
    'x',
]))

for i in reversed(range(ligand_kds.shape[0])):
    line_marker=plot_marker_styles[i]
    mark_every_n_points=NUM_POINTS_ON_XAXIS//10
    if line_marker=='x':
        mark_every_n_points=range((NUM_POINTS_ON_XAXIS//10)//2,NUM_POINTS_ON_XAXIS,NUM_POINTS_ON_XAXIS//10)
    if line_marker=='P':
        mark_every_n_points=range((NUM_POINTS_ON_XAXIS//10)//4,NUM_POINTS_ON_XAXIS,NUM_POINTS_ON_XAXIS//10)
    ax.plot(x_axis, y[i], 'k', label=plot_line_labels[i],
            marker=line_marker, markevery=mark_every_n_points, linewidth=1, markersize=7)
ax.set_xlabel(r"Inhibitor pK$_\mathrm{D}$")
ax.set_ylabel("Fraction ligand bound")
ax.legend()
ax.grid()
ax.title.set_text(r"Fraction ligand bound over a range of inhibitor K$_\mathrm{D}$s, [L$_0$]=100 nM, [I$_0$]=10 " +
                  r"$\mathrm{\mu}$M"+f"\n[P$_0$] = 8 nM.")
ax.set_xlim(3, 12)
ax.set_ylim(0, np.max(y)*1.1)
plt.show()
