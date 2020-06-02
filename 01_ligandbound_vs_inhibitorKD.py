import sys
from high_accuracy_binding_equations import *
from matplotlib import pyplot as plt
import numpy as np
 
 
NUM_POINTS_ON_X_AXIS=200
x_axis=np.linspace(3,9,NUM_POINTS_ON_X_AXIS)
inhibitor_kd_range=10**(-x_axis)*1e6
ligand_kds=np.array([0.001,0.5, 1,2,5,10])
protein_concs=np.array(calc_amount_a_i_from_kd(0.7, ligand_kds))
print(protein_concs)
ligand_kd_labels=[f"KDPL=1 nM, [P]={round(protein_concs[0],3)}", f"KDPL=500 nM, [P]={round(protein_concs[1],3)}", r"KDPL=1 $\mathrm{\mu}$M"+f", [P]={round(protein_concs[2],3)}",r"KDPL=2 $\mathrm{\mu}$M"+f", [P]={round(protein_concs[3],3)}",r"KDPL=5 $\mathrm{\mu}$M"+f", [P]={round(protein_concs[4],3)}",r"KDPL=10 $\mathrm{\mu}$M"+f", [P]={round(protein_concs[5],3)}",]
y=np.full((ligand_kds.shape[0],x_axis.shape[0]),np.nan)

for it_ligand_kds, ligand_kd in enumerate(ligand_kds):
    for it_inhibitor_kd_range, inhibitor_kd in enumerate(inhibitor_kd_range):
        lig_conc=ligand_kd*2
        i_conc=10
        p=protein_concs[it_ligand_kds]
        y[it_ligand_kds][it_inhibitor_kd_range]=competition_pl(**{'p':p,'l':lig_conc,'i':i_conc, 'kdpl':ligand_kd, 'kdpi':inhibitor_kd})/lig_conc
     
 
fig, ax = plt.subplots()
ax.set_xticklabels(["3 (mM)","4","5", r"6 ($\mathrm{\mu}$M)","7", "8", "9 (nM)"])
ax.set_xticks([3,4,5,6,7,8,9])
for i in reversed(range(ligand_kds.shape[0])): ax.plot(x_axis,y[i], label=ligand_kd_labels[i])
ax.set_xlabel("Inhibitor pKd")
ax.set_ylabel("Fraction ligand bound")
ax.legend()
ax.title.set_text("Fraction ligand bound over a range of inhibitor Kds, [L]=2xKDPL, [I]=10 "+r"$\mathrm{\mu}$M")
ax.set_xlim(3,9)
ax.set_ylim(0,0.75)
 
plt.show()