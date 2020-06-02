import sys
from high_accuracy_binding_equations import *
from matplotlib import pyplot as plt
import numpy as np

ligand_conc=0.01
ligand_kd=0.1
 
inhibitor_kd=1
inhibitor_conc=10

protein_conc=calc_amount_p(0.7,ligand_conc, ligand_kd)

print(competition_pl(**{'p':protein_conc,'l':ligand_conc,'i':0, 'kdpl':ligand_kd, 'kdpi':inhibitor_kd})/ligand_conc)

print(f"{inhibitor_kd=}, {ligand_kd=}, {protein_conc=}, ",)
print(round(competition_pl(**{'p':protein_conc,'l':ligand_conc,'i':inhibitor_conc, 'kdpl':ligand_kd, 'kdpi':inhibitor_kd})/ligand_conc,3))