"""
Plot grid varying L0 and target fraction ligand bound.

"""

import sys
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from high_accuracy_binding_equations import *

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# Parameters dictating range of simulation
DATAFILENAME="09-data"
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 1000 # Publication used 2000 pts along X
LIGAND_CONCS = [0.001, 0.010, 0.100]
TARGET_FLBS = [0.9,0.7,0.3]
INHIBITOR_CONC = 10
inhibitor_kds = np.array([0.001, 0.01, 0.1000001, 1, 10, 100])

if not Path(DATAFILENAME+".npz").exists():
    print("Data file doesnt exist, generating....")
    x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
    ligand_kd_range = 10**(-x_axis)*1e6  # We are working in µM, which is 1e-6.
    protein_concs = np.full((len(TARGET_FLBS), len(LIGAND_CONCS), len(ligand_kd_range)), np.nan)
    print("Reticulating splines...")
    for i_tflb, target_fraction_ligand_bound in enumerate(TARGET_FLBS):
        for i in range(protein_concs.shape[1]):
            protein_concs[i_tflb,i]=calc_amount_p(target_fraction_ligand_bound, LIGAND_CONCS[i], ligand_kd_range)
    print("Completed [P]s")
    y = np.full((len(TARGET_FLBS),len(LIGAND_CONCS), inhibitor_kds.shape[0], x_axis.shape[0]), np.nan)

    for i_tflb, target_fraction_ligand_bound in enumerate(TARGET_FLBS):
        for ligand_conc_i in range(len(LIGAND_CONCS)):
            for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
                for it_ligand_kd_range, ligand_kd in enumerate(ligand_kd_range):
                    p = protein_concs[i_tflb,ligand_conc_i,it_ligand_kd_range]
                    y[i_tflb,ligand_conc_i,it_inhibitor_kds,it_ligand_kd_range] = competition_pl(
                        **{'p': p, 'l': LIGAND_CONCS[ligand_conc_i], 'i': INHIBITOR_CONC, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/LIGAND_CONCS[ligand_conc_i]
            print(inhibitor_kd, end=",")
            sys.stdout.flush()
        print()


        print(f"Ligand conc {LIGAND_CONCS[ligand_conc_i]} done")
    np.savez_compressed(DATAFILENAME, x_axis=x_axis, ligand_kd_range=ligand_kd_range,protein_concs=protein_concs,y=y)
    print("Made datafile")

print("Loading datafile")
loaded_data=np.load(DATAFILENAME+".npz")
x_axis=loaded_data['x_axis']
ligand_kd_range=loaded_data['ligand_kd_range']
protein_concs=loaded_data['protein_concs']
y=loaded_data['y']
#
fig, ax = plt.subplots(len(LIGAND_CONCS), len(TARGET_FLBS), figsize=(6,8), sharex='col', sharey='row')
ax[-1,-1].set_xticklabels(["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax[-1,-1].set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

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
for i_tflb, tflb in enumerate(TARGET_FLBS):
    for i_ligand_concs in range(len(LIGAND_CONCS)):
        for i in reversed(range(inhibitor_kds.shape[0])):
            ax[i_tflb, i_ligand_concs].plot(x_axis, y[i_tflb,i_ligand_concs,i], 'k', label=plot_line_labels[i],
                    marker=plot_marker_styles[i], markevery=NUM_POINTS_ON_XAXIS//5, linewidth=1)
        
ax[0,0].legend()
fig.suptitle("Protein-ligand signal over a range of ligand K$_\mathrm{D}$s, [I$_0$]=10 $\mathrm{\mu}$M,\n")
# For publication supporting information figures, change layout to the following:
# Left = 0.09, Bottom = 0.05, Right = 0.96, Top = 0.94, Wspace = 0.2, Hspace = 0.16
for i in range(len(TARGET_FLBS)):
    for j in range(len(LIGAND_CONCS)):
         ax[i,j].set_xlim(3, 12)
         #ax[i,j].grid()
plt.tight_layout()
plt.show()
