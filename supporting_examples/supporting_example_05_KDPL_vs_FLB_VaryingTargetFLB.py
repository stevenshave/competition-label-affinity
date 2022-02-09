"""
Plot competition experiment sensitivity, varying the target fraction bound.

Generates a plot used in the Supporting Information section of the manuscript
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
NUM_POINTS_ON_XAXIS = 40 # Publication used 2000 pts along X
TARGET_FRACTION_L_BOUND_LIST = [0.7,0.5,0.3,0.1]
LIGAND_CONC = 0.01
INHIBITOR_CONC = 10

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10**(-x_axis)*1e6  # We are working in µM, which is 1e-6.
inhibitor_kds = np.array([0.001, 0.01, 0.1000001, 1, 10, 100])
protein_concs = np.full((len(TARGET_FRACTION_L_BOUND_LIST), len(ligand_kd_range)), np.nan)
print("Reticulating splines...")

for i in range(protein_concs.shape[0]):
    protein_concs[i]=calc_amount_p(TARGET_FRACTION_L_BOUND_LIST[i], 0.01, ligand_kd_range)
print("Completed [P]s")

y = np.full((len(TARGET_FRACTION_L_BOUND_LIST), inhibitor_kds.shape[0], x_axis.shape[0]), np.nan)

for tfb_index in range(len(TARGET_FRACTION_L_BOUND_LIST)):
    for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
        for it_ligand_kd_range, ligand_kd in enumerate(ligand_kd_range):
            p = protein_concs[tfb_index,it_ligand_kd_range]
            y[tfb_index,it_inhibitor_kds,it_ligand_kd_range] = competition_pl(
                **{'p': p, 'l': LIGAND_CONC, 'i': INHIBITOR_CONC, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/LIGAND_CONC
    print(f"Target fraction ligand bound {TARGET_FRACTION_L_BOUND_LIST[tfb_index]} done")
    

fig, ax = plt.subplots(4,1, figsize=(8, 6*len(TARGET_FRACTION_L_BOUND_LIST)))
ax[-1].set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax[-1].set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

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

for plotindex in range(len(TARGET_FRACTION_L_BOUND_LIST)):
    for i in reversed(range(inhibitor_kds.shape[0])):
        ax[plotindex].plot(x_axis, y[plotindex,i], 'k', label=plot_line_labels[i],
                marker=plot_marker_styles[i], markevery=NUM_POINTS_ON_XAXIS//20, linewidth=1)
        ax[plotindex].set_xlim(3, 12)
        ax[plotindex].set_ylim(0, TARGET_FRACTION_L_BOUND_LIST[plotindex]*1.1)
        ax[plotindex].set_ylabel("Fraction ligand bound")
        ax[plotindex].title.set_text(f"Target fraction ligand bound = {TARGET_FRACTION_L_BOUND_LIST[plotindex]}")
ax[0].legend()
ax[-1].set_xlabel(r"Ligand pK$_\mathrm{D}$")
fig.suptitle("Protein-ligand signal over a range of ligand K$_\mathrm{D}$s, [L$_0$]=10 nM, [I$_0$]=10 " +
                  r"$\mathrm{\mu}$M")

# For publication supporting information figures, change layout to the following:
# Left = 0.09, Bottom = 0.05, Right = 0.96, Top = 0.94, Wspace = 0.2, Hspace = 0.16

plt.show()
