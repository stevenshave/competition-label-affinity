
# competition-label-affinity 
[![DOI](https://zenodo.org/badge/264960530.svg)](https://zenodo.org/badge/latestdoi/264960530)

A tool for investigating the impact of labelled ligand affinity in competition experiments.  Code accompanies the paper "Optimum ligand affinity for competition-based primary screens" by Shave et. al.

Can be installed from pypi via
``` pip install claffinity ```

![LigandpKD vs FLB](https://github.com/stevenshave/competition-label-affinity/blob/7c3a184d4cc8127fcf22e1f125721a1b66ba174f/Figure01-0.7.png "LigandpKD vs FLB")


# Example programs:
A selection of programs are available in this archive prefixed with 'example_'

A breakdown of their use in simulating and understanding the imact of ligand affinity in competition assays follows:

### example_01_plot_protein_needed_vs_ligand_kd.py
 - Plot ligand pKD vs conc protein needed to achieve desired (default = .7) fraction ligand bound.
### example_01_plot_singal_vs_ligand_kd_fixed_inhibitor.py
 - Plot ligand pKD vs fraction ligand bound
### example_01_simulate_competition_experiment.py
 - Simultate a simple single point from a competition experiment.

### example_02_inhibitorKD_vs_fractionBound.py
 - Simulate inhibitor pKD vs fraction ligand bound
### example_02_ligandKD_vs_fractionBound.py
 - Plot ligand pKD vs fraction ligand bound in the presence of a fixed concentration of inhibitor with fixed KD.
 ### example_03_ligandKD_vs_fractionBoundAsPercentSignalReduction.py
- Plot ligand pKD vs percent signal reduction in the presence of a fixed concentration of inhibitor with fixed KD.

### example_04_missed_inhibitors.py
- Simulate real world example whereby inhibitors would be missed in a primary screen using high affinity ligands.

## Supporting example programs
Some additional example application of the simulation techniques outlined in the paper are shown below, including code used in supporting information figure generation, the generation of animations and the SI matterial video.
### supporting_example_02_inhibitorKD_vs_fractionBound_animation.py
- Generate animation of inhibitor pKD vs fraction ligand bound over a range of ligand pKDs in the animation.
### supporting_example_02_inhibitorKD_vs_fractionBound_interactive.py
- Interactive plot of inhibitor pKD vs fraction ligand bound over a range of changable ligand pKDs.

### supporting_example_02_ligandKD_vs_fractionBound_animation.py
 - Generate animation for ligand pKD vs fraction ligand bound, varying ligand KD over time.

### supporting_example_02_ligandKD_vs_fractionBound_interactive.py
 - Interactive plot of ligand pKD vs fraction ligand bound, varying inhibitor KD.

### supporting_example_05_KDPL_vs_FLB_VaryingTargetFLB.py
- Plot pKD of ligand vs fraction ligand bound, varying the target fraction ligand bound from the default 0.7.

### supporting_example_06_KDPL_vs_FLB_VaryingL0.py
- Plot pKD of ligand vs fraction ligand bound, varying the amount of ligand present.

### supporting_example_07_KDPL_vs_FLB_VaryingI0.py
- Plot pKD of ligand vs fraction ligand bound, varying the amount of inhibitor present.

### supporting_example_08_HuangPlot.py
- Reproduce the Huang plot *(Huang, X., Fluorescence polarization competition assay: the range of resolvable inhibitor potency is limited by the affinity of the fluorescent ligand. Journal of biomolecular screening 2003, 8 (1), 34-38.)*


## Requirements
Code developed using python 3.7.1 but should work with any Python version 3.6 or greater. The following packages are also required
- matplotlib
- numpy>=1.15
- pandas>=1.2.2
- mpmath>=1.1.0
- progressbar2

