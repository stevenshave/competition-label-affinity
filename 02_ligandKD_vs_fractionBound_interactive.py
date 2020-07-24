"""
Produce plot of competition experiment sensitivity

Generates a plot used in the manuscript "Revisiting labelled ligand affinity
in competition experiments" by Shave et.al.
"""

import sys
from matplotlib import pyplot as plt
from pathlib import Path
plt.rcParams['animation.ffmpeg_path'] = Path("C:\\Users\\steve\\Downloads\\ffmpeg-20200716-d11cc74-win64-static\\bin\\ffmpeg.exe")
from matplotlib import animation
from matplotlib.widgets import Slider
import sys
import numpy as np
from high_accuracy_binding_equations import *



# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
TARGET_FRACTION_L_BOUND = 0.7
NUM_INHIBITOR_KDS = 35
NUM_LIGAND_KDS = 200

x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END,NUM_LIGAND_KDS)
inhibitor_kds = 10**(-np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_INHIBITOR_KDS))  
ligand_kds =  10**(-np.linspace(XAXIS_BEGINNING,XAXIS_END, num=NUM_LIGAND_KDS))

protein_concs = np.array(calc_amount_p(
    TARGET_FRACTION_L_BOUND, 10e-9, ligand_kds))
DATAFILENAME="02-data.npz"

y = np.full((ligand_kds.shape[0], inhibitor_kds.shape[0]), np.nan)
if not Path(DATAFILENAME).exists():
    print("Regenerating data file")
    for it_ligand_kds, ligand_kd in enumerate(ligand_kds):
        print(f"{it_ligand_kds}/{ligand_kds.shape[0]}")
        
        for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
            lig_conc = 10e-9
            i_conc = 10e-6
            p = protein_concs[it_ligand_kds]
            y[it_ligand_kds][it_inhibitor_kds] = competition_pl(
                **{'p': p, 'l': lig_conc, 'i': i_conc, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/lig_conc
    np.savez_compressed(DATAFILENAME, xaxis_beginning=XAXIS_BEGINNING, xaxis_end=XAXIS_END, ligand_kds=ligand_kds,protein_concs=protein_concs,y=y, inhibitor_kds=inhibitor_kds)

print("Loading datafile")
loaded_data=np.load(DATAFILENAME)
XAXIS_BEGINNING=loaded_data['xaxis_beginning']
XAXIS_END=loaded_data['xaxis_end']
ligand_kds=loaded_data['ligand_kds']
inhibitor_kds=loaded_data['inhibitor_kds']
protein_concs=loaded_data['protein_concs']
x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END, num=ligand_kds.shape[0])
y=loaded_data['y']


fig, ax = plt.subplots(2,1, figsize=(7.204724, 5.09424929292), gridspec_kw={'height_ratios':[10,1]})
line,=ax[0].plot(x_axis, y[:,0], 'k',  linewidth=1)

pkd_to_animation_frame=lambda x: int((x-XAXIS_BEGINNING)/(XAXIS_END-XAXIS_BEGINNING)*inhibitor_kds.shape[0])
def update(num):
    intermediate=pkd_to_animation_frame(num)
    if intermediate==inhibitor_kds.shape[0]: intermediate=inhibitor_kds.shape[0]-1
    line.set_data(x_axis, y[:,int(intermediate)])
    fig.canvas.draw_idle()

update(6)



ligand_kd_slider = Slider(ax[1], r"Inhibitor pK$_\mathrm{D}$", XAXIS_BEGINNING, XAXIS_END, valinit=6)
ligand_kd_slider.on_changed(update)






ax[0].set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax[0].set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))






ax[0].set_xlabel(r"Ligand pK$_\mathrm{D}$")
ax[0].set_ylabel("Fraction ligand bound")
ax[0].grid()
ax[0].title.set_text(r"Protein-ligand signal over a range of ligand K$_\mathrm{D}$s, [L]=10 nM, [I]=10 " +
                  r"$\mathrm{\mu}$M"+f"\nTarget fraction ligand bound without inhibitor = {TARGET_FRACTION_L_BOUND}")
ax[0].set_xlim(3, 12)
ax[0].set_ylim(0, TARGET_FRACTION_L_BOUND*1.1)
plt.tight_layout()
plt.show()
