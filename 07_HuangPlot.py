"""
Produce plot similar to Huang

Generates a plot used in the supporting information of the manuscript
"Revisiting labelled ligand affinity in competition experiments" by Shave
et.al.
"""

import sys
from matplotlib import pyplot as plt
import numpy as np
from high_accuracy_binding_equations import *

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume nM for all
# concentrations bellow.

# Parameters dictating range of simulation
XAXIS_BEGINNING = 0  # 10**0 ==1 nM
XAXIS_END = 5  # 10**5 = 100 µM
NUM_POINTS_ON_XAXIS = 2000  # Publication used 2000 pts along X
TARGET_FRACTION_L_BOUND = 0.7
LIGAND_KDs=[10,100,1000,10000]
LIGAND_CONC=10

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
inhibitor_kds = 10**(x_axis)  # We are working in µM, which is 1e-6.

protein_conc_for_ligand_kd = np.full((len(LIGAND_KDs)),np.nan)
for i,ligand_kd in enumerate(LIGAND_KDs):
    protein_conc_for_ligand_kd[i]=calc_amount_p(TARGET_FRACTION_L_BOUND, LIGAND_CONC, ligand_kd)
print(protein_conc_for_ligand_kd)
y = np.full((protein_conc_for_ligand_kd.shape[0], inhibitor_kds.shape[0]), np.nan)

for it_protein_conc, protein_conc in enumerate(protein_conc_for_ligand_kd):
    for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
        y[it_protein_conc][it_inhibitor_kds]=calc_i_for_fractionl_bound(protein_conc,LIGAND_CONC,LIGAND_KDs[it_protein_conc],inhibitor_kd,0.7/2)

plot_line_labels = [
    r'K$_\mathrm{D}$PL=10 nM',
    r'K$_\mathrm{D}$PL=100 nM',
    r'K$_\mathrm{D}$PL=1 µM',
    r'K$_\mathrm{D}$PL=10 µM',
]
plot_marker_styles = [
    's',
    '^',
    'v',
    'D',
]
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xticklabels(
    ["1 nM", "10 nM", "100 nM", "1 µM", "10 µM", "100 µM"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

for i in reversed(range(y.shape[0])):
    ax.plot(x_axis, y[i], 'k', marker=plot_marker_styles[i], markevery=NUM_POINTS_ON_XAXIS//10, linewidth=1,label=plot_line_labels[i])
ax.set_xlabel(r"K$_\mathrm{I}$ (KDPI)")
ax.set_ylabel(r"IC$_{50}$ ([I]) (nM) ")
ax.legend()
ax.grid()
plt.yscale('log')
ax.title.set_text(r"IC$_{50}$ versus K$_i$, L$_0$ = 10 nM, target fraction ligand bound = "+str(TARGET_FRACTION_L_BOUND)) 
ax.set_xlim(XAXIS_BEGINNING, XAXIS_END)
plt.show()
