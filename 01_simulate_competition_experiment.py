"""Simulation 1:1:1 comptition binding"""

import numpy as np
from high_accuracy_binding_equations import *

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# First, lets simulate a few single points from three different systems
# p, l and i are protrin, ligand and inhibitor concentrations respectively
# kdpl is the dissociation constant (KD) of the protein-ligand interaction
# kdpi is the dissociation constant (KD) of the protein-inhibitor interaction
# We can either expand the dictionary with ** as shown in the example with
# system1, or we can pass arguments to competition_pl with the following
# singature: competition_pl(p, l , i, kdpl, kdpi)

system1={"p":1, "l":2, "i":10, "kdpl":0.1, "kdpi":10}
pl_conc=competition_pl(**system1)
print(f"pl_conc = {round(pl_conc,4)}, fraction ligand bound = {round(pl_conc/system1['l'],4)}")






