import sys
from matplotlib import pyplot as plt
import matplotlib.ticker as plticker
import numpy as np
from high_accuracy_binding_equations import *

x_axis=np.linspace(3,9,400)
print(x_axis)
ligand_kd_range=10**(-x_axis)*1e6
inhibitor_kds=np.array([0.001,0.01,0.1000001,1,10,100])
protein_concs=np.array(calc_amount_p(0.7, 0.01, ligand_kd_range))
y=np.full((inhibitor_kds.shape[0],x_axis.shape[0]),np.nan)
 
for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
    for it_ligand_kd_range, ligand_kd in enumerate(ligand_kd_range):
        lig_conc=0.01
        i_conc=10
        p=protein_concs[it_ligand_kd_range]
        y[it_inhibitor_kds][it_ligand_kd_range]=competition_pl(**{'p':p,'l':lig_conc,'i':i_conc, 'kdpl':ligand_kd, 'kdpi':inhibitor_kd})/lig_conc
     
 
fig, ax = plt.subplots()
ax.set_xticklabels(["3 (mM)","4","5", r"6 ($\mathrm{\mu}$M)","7", "8", "9 (nM)"])
ax.set_xticks([3,4,5,6,7,8,9])
for i in reversed(range(inhibitor_kds.shape[0])): ax.plot(x_axis,y[i], label=f"KDPI = {inhibitor_kds[i]:4.4}"+r" $\mathrm{\mu}$M")
#   ax.plot(x_axis, protein_concs/10, label="[P]")
ax.set_xlabel("Ligand pKd")
ax.set_ylabel("Fraction ligand bound")
ax.legend()
ax.title.set_text("Fraction ligand bound over a range of ligand Kds, [L]=10 nM, [I]=10 "+r"$\mathrm{\mu}$M")
ax.set_xlim(3,12)
ax.set_ylim(0,0.75)
 
plt.show()