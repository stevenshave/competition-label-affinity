"""
Produce plot of competition experiment sensitivity

Generates a plot used in the manuscript "Identification of optimum ligand affinity
for competition-based primary screens" by Shave et.al.
"""

import sys
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
from pathlib import Path
plt.rcParams['animation.ffmpeg_path'] = Path("C:\\Program Files\\ffmpeg\\bin\\ffmpeg.exe")
from matplotlib import animation
import sys
import numpy as np
from claffinity.high_accuracy_binding_equations import *
from matplotlib import cm


# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
TARGET_FRACTION_L_BOUND = 0.7
NUM_INHIBITOR_KDS = 1000
NUM_LIGAND_KDS = 1000

x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END,NUM_INHIBITOR_KDS)
inhibitor_kds = 10**(-np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_INHIBITOR_KDS))  
ligand_kds =  10**(-np.linspace(XAXIS_BEGINNING,XAXIS_END, num=NUM_LIGAND_KDS))

protein_concs = np.array(calc_amount_p(
    TARGET_FRACTION_L_BOUND, 10e-9, ligand_kds))
DATAFILENAME="02-data.npz"

flb = np.full((ligand_kds.shape[0], inhibitor_kds.shape[0]), np.nan)
if not Path(DATAFILENAME).exists():
    print("Regenerating data file")
    for it_ligand_kds, ligand_kd in enumerate(ligand_kds):
        print(f"{it_ligand_kds}/{ligand_kds.shape[0]}")
        
        for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
            lig_conc = 10e-9
            i_conc = 10e-6
            p = protein_concs[it_ligand_kds]
            flb[it_ligand_kds][it_inhibitor_kds] = competition_pl(
                **{'p': p, 'l': lig_conc, 'i': i_conc, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/lig_conc
    np.savez_compressed(DATAFILENAME, xaxis_beginning=XAXIS_BEGINNING, xaxis_end=XAXIS_END, ligand_kds=ligand_kds,protein_concs=protein_concs,y=flb, inhibitor_kds=inhibitor_kds)

print("Loading datafile")
loaded_data=np.load(DATAFILENAME)
XAXIS_BEGINNING=loaded_data['xaxis_beginning']
XAXIS_END=loaded_data['xaxis_end']
ligand_kds=loaded_data['ligand_kds']
inhibitor_kds=loaded_data['inhibitor_kds']
protein_concs=loaded_data['protein_concs']
x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END, num=ligand_kds.shape[0])
y_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END, num=inhibitor_kds.shape[0])
flb=loaded_data['y']

#fig, ax = plt.subplots(2,1, figsize=(7.204724, 5.09424929292), gridspec_kw={'height_ratios':[10,1]})
fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(7.204724, 5.09424929292))

X, Y = np.meshgrid(x_axis, y_axis)
ax.plot_surface(X, Y, flb, linewidth=1, cmap=cm.coolwarm,antialiased=False, rstride=10,cstride=10, alpha=0.8)



ax.set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

ax.set_yticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax.set_yticks(range(XAXIS_BEGINNING, XAXIS_END+1))

find_index_of_closest_to=lambda val, nparray: np.abs(nparray-val).argmin()
optimal_y=np.full_like(y_axis, y_axis[find_index_of_closest_to(6.975, y_axis)])
optimal_flb=flb[find_index_of_closest_to(6.975, x_axis), :]
ax.scatter(x_axis, optimal_y, optimal_flb, 'k', s=10, alpha=1)


ax.set_ylabel(r"Ligand pK$_\mathrm{D}$")
ax.set_xlabel("Inhibitor pK$_\mathrm{D}$")
ax.set_zlabel("Fraction ligand bound")
ax.grid()
ax.title.set_text(r"Fraction ligand bound over a range of ligand and inhibitor K$_\mathrm{D}$s, [L$_0$]=10 nM, [I$_0$]=10 " +
                  r"$\mathrm{\mu}$M"+f"\nTarget fraction ligand bound without inhibitor = {TARGET_FRACTION_L_BOUND}")
#ax.set_xlim(3, 12)
#ax.set_ylim(0, TARGET_FRACTION_L_BOUND*1.1)
plt.show()
