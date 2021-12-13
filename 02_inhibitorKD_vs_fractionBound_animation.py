"""
Produce plot of competition experiment sensitivity

Generates a plot used in the manuscript "Identification of optimum ligand affinity
for competition-based primary screens" by Shave et.al.
"""

import sys
from matplotlib import pyplot as plt
from pathlib import Path
plt.rcParams['animation.ffmpeg_path'] = Path("C:\\Users\\steve\\Downloads\\ffmpeg-20200716-d11cc74-win64-static\\bin\\ffmpeg.exe")
from matplotlib import animation
import sys
import numpy as np
from high_accuracy_binding_equations import *



# Parameters dictating range of simulation
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
TARGET_FRACTION_L_BOUND = 0.7
NUM_INHIBITOR_KDS = 1000
NUM_LIGAND_KDS = 1000
USE_VIDEO_SIZING=True


figure_size=(7.204724, 5.09424929292)
plot_title_size=16
axis_label_size=14
tick_label_font_size=12
if USE_VIDEO_SIZING:
    figure_size=(19.2,10.8)
    plot_title_size=28
    axis_label_size=20
    tick_label_font_size=20

x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END,NUM_INHIBITOR_KDS)
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
x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END, num=inhibitor_kds.shape[0])
y=loaded_data['y']


#fig.set_size_inches(*figure_size, forward = False)
fig, ax = plt.subplots(2,1, figsize=figure_size, gridspec_kw={'height_ratios':[10,1]}, sharex=True)
line,=ax[0].plot(x_axis, y[0], 'k',  linewidth=1)


for tick in ax[1].xaxis.get_major_ticks():
    tick.label.set_fontsize(tick_label_font_size) 
for tick in ax[1].yaxis.get_major_ticks():
    tick.label.set_fontsize(tick_label_font_size) 
for tick in ax[0].yaxis.get_major_ticks():
    tick.label.set_fontsize(tick_label_font_size) 


def update(num, x, y, line):
    line.set_data(x_axis, y[int(num)])
    for bar in ax[1].containers:
        bar.remove()
    ax[1].barh(y=r"Ligand pK$_\mathrm{D}$", width=np.linspace(XAXIS_BEGINNING,XAXIS_END, num=ligand_kds.shape[0])[num], color='k')
    fig.canvas.draw_idle()
    return line

ani = animation.FuncAnimation(fig, update, range(len(ligand_kds))[::4], fargs=[x_axis, y, line],interval=5, blit=False, repeat_delay=2000,)

ax[0].set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)", "10", "11", "12 (pM)"])
ax[0].set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))

ax[0].set_xlabel(r"Inhibitor pK$_\mathrm{D}$", fontsize=axis_label_size)
ax[0].set_ylabel("Fraction ligand bound", fontsize=axis_label_size)
ax[0].grid()
fig.suptitle(r"Protein-ligand signal over a range of inhibitor K$_\mathrm{D}$s, [L$_0$]=10 nM, [I$_0$]=10 " +
                  r"$\mathrm{\mu}$M"+f"\nTarget fraction ligand bound without inhibitor = {TARGET_FRACTION_L_BOUND}",fontsize=plot_title_size)
ax[0].set_xlim(3, 12)
ax[0].set_ylim(0, TARGET_FRACTION_L_BOUND*1.1)
fig.tight_layout(rect=[0.04,0,1,0.9])

ani.save('animation_inhib_vs_fb_fast.mp4', writer = "ffmpeg", extra_args=['-vcodec', 'libx264'],fps=30)
plt.show()
