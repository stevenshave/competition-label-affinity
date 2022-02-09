


import numpy as np
import random
import numpy as np
import pandas as pd
import sys

from claffinity.high_accuracy_binding_equations import *
XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 12  # pKD of 12 is pM
NUM_LIGAND_KDS = 200


x_axis=np.linspace(XAXIS_BEGINNING, XAXIS_END,NUM_LIGAND_KDS)
ligand_kds =  10**(-np.linspace(XAXIS_BEGINNING,XAXIS_END, num=NUM_LIGAND_KDS))

y = np.full((ligand_kds.shape[0]), np.nan)

TFLBs=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
pl0s=[1e-12]
while pl0s[-1]!=1e-3:
    pl0s.append(pl0s[-1]*5)
    pl0s.append(pl0s[-1]*2)
pl0s=-np.log10(pl0s[::-1])


pi0s=pl0s[0:9]
n=0
for TFLB in TFLBs:
    frame=pd.DataFrame(index=pi0s,columns=pl0s)
    for pl0 in pl0s:
        for pi0 in pi0s:
            ligand_conc=10**-pl0#10**-np.random.normal(5,0.5)
            inhibitor_conc=10**-pi0#10**-np.random.normal(5,0.5)
            inhibitor_kd=10**-6
            protein_conc=np.array(calc_amount_p(TFLB, ligand_conc, ligand_kds))
            # Simulate
            for it_ligand_kd_range, ligand_kd in enumerate(ligand_kds):
                y[it_ligand_kd_range] = competition_pl(**{'p': protein_conc[it_ligand_kd_range], 'l': ligand_conc, 'i': inhibitor_conc, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/ligand_conc
            frame.loc[pi0,pl0]=-np.log10(ligand_kds[np.argmin(y)])
            #frame = frame.append({"TFLB":TFLB, "p[L0]":pl0, "p[I0]":pi0, "pKDPLmin":-np.log10(ligand_kds[np.argmin(y)])}, ignore_index=True)
            n+=1
            print(frame)
            print(f"n={n}/{len(TFLBs)*len(pl0s)*len(pi0s)}")
    frame.to_csv(f"lookup_table{TFLB}.csv")
