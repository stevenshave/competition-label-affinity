from typing import Optional, Union, List, Tuple
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
import pandas as pd
import progressbar

from .high_accuracy_binding_equations import (
    calc_amount_p,
    calc_kdpi_for_fractionl_bound,
    calc_i_for_fractionl_bound,
    competition_pl,
)
from math import floor, ceil


class CompetitionLabelAffinity:
    

    def __init__(
        self,
        p: Optional[float] = None,
        l: Optional[float] = None,
        i: Optional[float] = None,
        kdpl: Optional[float] = None,
        kdpi: Optional[float] = None,
        x_axis_resolution: float = 40,
    ):
        self.pkd_label_lookup = {
            "1": "1 (M)",
            "3": "3 (mM)",
            "6": r"6 ($\mathrm{\mu}$M)",
            "9": "9 (nM)",
            "12": "12 (pM)",
            "15": "15 (fM)",
        }
        self.plot_marker_styles=list(reversed([
            '',
            '+',
            'x',
            "P",
            "s",
            "8",
            "*",
            6,
            7
        ]))

        self.p = p
        self.l = l
        self.i = i
        self.kdpl = kdpl
        self.kdpi = kdpi
        self.x_axis_resolution = x_axis_resolution
        self.kd_str=r'K$_\mathrm{D}$'

    def float_to_prettyprint_conc(self, f:float)->str:
        if f<=1e-9:
            return f"{f/1e-12:.2f} pM"
        if f<=1e-6:
            return f"{f/1e-9:.2f} nM"
        if f<=1e-3:
            return f"{f/1e-6:.2f} µM"
        return f"{f} M"

    def set_x_ticks_and_labels(self, axis, pKD_begin: float, pKD_end: float):
        if isinstance(pKD_begin, float):
            pKD_begin = int(floor(pKD_begin))
        if isinstance(pKD_end, float):
            pKD_end = int(ceil(pKD_end))

        axis.set_xticks(range(pKD_begin, pKD_end + 1))

        xtick_labels = [
            l if l not in self.pkd_label_lookup else self.pkd_label_lookup[l]
            for l in [f"{pkd}" for pkd in np.arange(pKD_begin, pKD_end + 1, 1)]
        ]
        axis.set_xticklabels(xtick_labels)
        axis.set_xlim(pKD_begin, pKD_end)

    def get_plot_line_labels(self, kds:Union[List[str], Tuple[str]]):
        line_labels=[]
        for kd in kds:
            if kd in self.pkid_label_lookup:
                line_labels.append(self.line_labels_lookup[str(kd)])
            else:
                line_labels.append(f"{kd}")
        return line_labels

    def plot_protein_needed_vs_ligand_kd(
        self,
        pKD_begin: float = 3,
        pKD_end: float = 12,
        l: float = 10e-9,
        target_fraction_ligand_bound=0.7,
    ):
        """
            Produce plot of amount of protein needed over a range of ligand KDs

        Generates a plot of protein needed to reach a desired fraction bound
        in a competition experiment, dependant on ligand KD.  For the manuscript
        "Identification of optimum ligand affinity for competition-based primary screens"
        by Shave et.al.

        Args:
            pKD_begin (float, optional): [description]. Defaults to 3.
            pKD_end (float, optional): [description]. Defaults to 12.
            ligand_conc (float, optional): [description]. Defaults to 10e-9.
            target_fraction_ligand_bound (float, optional): [description]. Defaults to 0.7.
            x_axis_resolution (int, optional): [description]. Defaults to 1000.
        """
        x_axis = np.linspace(pKD_begin, pKD_end, self.x_axis_resolution)
        ligand_kd_range = 10 ** (-x_axis)
        protein_concs = np.array(calc_amount_p(target_fraction_ligand_bound, l, ligand_kd_range))
        fig, ax = plt.subplots(figsize=(8, 6))
        self.set_x_ticks_and_labels(ax, pKD_begin, pKD_end)
        ax.plot(x_axis, protein_concs, "k")
        ax.set_xlabel(f"Ligand p{self.kd_str}")
        ax.set_ylabel(r"[P$_0$] (M)")
        plt.yscale("log")
        ax.grid()
        ax.title.set_text(
            f"Protein required for {target_fraction_ligand_bound} fraction ligand bound vs labeled ligand {self.kd_str}, [L$_0$]={self.float_to_prettyprint_conc(l)}"
        )
        ax.set_xlim(pKD_begin, pKD_end)
        plt.show()

    def plot_signal_vs_ligand_kd_fixed_i(
        self,
        pKD_begin: float = 3,
        pKD_end: float = 12,
        p: float = 10e-6,
        l: float = 10e-9,
        i: float = 10e-6,
        kdpi: float = 10e-6,
        target_fraction_ligand_bound=0.7,
    ):
        """Produce plot of competition experiment sensitivity with fixed lignand

        Generates a supporting plot of signal (PL) in a competition experiment vs
        ligand KD for a fixed concentration and KD of inhibitor "Identification
        of optimum ligand affinity for competition-based primary screens" by Shave et.al."""

        x_axis = np.linspace(pKD_begin, pKD_end, self.x_axis_resolution)
        ligand_kd_range = 10 ** (-x_axis)
        y = np.full((x_axis.shape[0]), np.nan)
        for idx in range(ligand_kd_range.shape[0]):
            y[idx] = competition_pl(**{"p": p, "l": l, "i": i, "kdpl": ligand_kd_range[idx], "kdpi": kdpi}) / l

        fig, ax = plt.subplots(figsize=(8, 6))

        ax.plot(x_axis, y, "k")
        self.set_x_ticks_and_labels(ax, pKD_begin, pKD_end)
        ax.set_xlabel(f"Ligand p{self.kd_str}")
        ax.set_ylabel("Fraction ligand bound")
        ax.title.set_text(
            f"Competition experiment fraction ligand bound over a range of ligand {self.kd_str}s"
            + "\n"
            + f"[P]={self.float_to_prettyprint_conc(p)}, [L]={self.float_to_prettyprint_conc(l)}, [I]={self.float_to_prettyprint_conc(i)}, inhibitor {self.kd_str}={self.float_to_prettyprint_conc(kdpi)}"
        )
        ax.grid()
        ax.set_ylim(0, 1.05)
        plt.show()
    def calc_amount_p(self,fraction_bound, l, kdax):
        return float(calc_amount_p(fraction_bound,l,kdax).real)

    def single_point_competition_readout(self, p: float, l: float, i: float, kdpl: float, kdpi: float):
        return float(competition_pl(p, l, i, kdpl, kdpi).real)

    def plot_inhibitor_KD_vs_FLB(self,
        pKD_begin: float = 3,
        pKD_end: float = 12,
        l: float = 10e-9,
        i: float = 10e-6,
        kdpl: Union[float, List[float], Tuple[float]] = [1e-12,10e-12,100e-12,1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6],
        target_fraction_ligand_bound=0.7,
    ):
        if isinstance(kdpl, (float, int)):
            kdpl=[kdpl]
        # Parameters dictating range of simulation
        x_axis = np.linspace(pKD_begin, pKD_end, self.x_axis_resolution)
        inhibitor_kd_range = 10**(-x_axis)
        if isinstance(kdpl, float):
            kdpl=[kdpl]
        kdpl = np.array(kdpl)
        protein_concs = np.array(calc_amount_p(
            target_fraction_ligand_bound, l, kdpl))
        y = np.full((kdpl.shape[0], x_axis.shape[0]), np.nan)

        for it_ligand_kd, ligand_kd in enumerate(kdpl):
            print(f"Generating: {it_ligand_kd+1}/{len(kdpl)}")
            for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kd_range):
                p = protein_concs[it_ligand_kd]
                y[it_ligand_kd][it_inhibitor_kds] = competition_pl(
                    **{'p': p, 'l': l, 'i': i, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/l

        fig, ax = plt.subplots(figsize=(8, 6))
        self.set_x_ticks_and_labels(ax, pKD_begin, pKD_end)
        
        for idx in reversed(range(kdpl.shape[0])):
            line_marker=self.plot_marker_styles[idx]
            mark_every_n_points=self.x_axis_resolution//10
            if line_marker=='x':
                mark_every_n_points=range((self.x_axis_resolution//10)//2,self.x_axis_resolution,self.x_axis_resolution//10)
            if line_marker=='P':
                mark_every_n_points=range((self.x_axis_resolution//10)//4,self.x_axis_resolution,self.x_axis_resolution//10)
            ax.plot(x_axis, y[idx], 'k', label=self.float_to_prettyprint_conc(kdpl[idx]),
                    marker=line_marker, markevery=mark_every_n_points, linewidth=1, markersize=7)
        ax.set_xlabel(f"Inhibitor p{self.kd_str}")
        ax.set_ylabel("Fraction ligand bound")
        ax.legend()
        ax.grid()
        ax.title.set_text(f"Fraction ligand bound over a range of inhibitor {self.kd_str}s, [L$_0$]={self.float_to_prettyprint_conc(l)}, [I$_0$]={self.float_to_prettyprint_conc(i)}"
            +f"\nTarget fraction ligand bound without inhibitor = {target_fraction_ligand_bound}")
        ax.set_xlim(pKD_begin, pKD_end)
        ax.set_ylim(0, target_fraction_ligand_bound*1.1)
        plt.show()

    def plot_ligand_KD_vs_FLB(self,
        pKD_begin: float = 3,
        pKD_end: float = 12,
        l: float = 10e-9,
        i: float = 10e-6,
        kdpi: Union[float, List[float], Tuple[float]] = [1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6],
        target_fraction_ligand_bound=0.7,
    ):
        if isinstance(kdpi, (float, int)):
            kdpi=[kdpi]
        x_axis = np.linspace(pKD_begin, pKD_end, self.x_axis_resolution)
        ligand_kd_range = 10**(-x_axis)  
        inhibitor_kds = np.array(kdpi)
        protein_concs = np.array(calc_amount_p(
            target_fraction_ligand_bound, l, ligand_kd_range))
        y = np.full((inhibitor_kds.shape[0], x_axis.shape[0]), np.nan)

        for it_inhibitor_kds, inhibitor_kd in enumerate(inhibitor_kds):
            print("Generating : "+str(it_inhibitor_kds+1)+" / "+str(len(inhibitor_kds)))
            for it_ligand_kd_range, ligand_kd in enumerate(ligand_kd_range):
                p = protein_concs[it_ligand_kd_range]
                y[it_inhibitor_kds][it_ligand_kd_range] = competition_pl(
                    **{'p': p, 'l': l, 'i': i, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/l

        fig, ax = plt.subplots(figsize=(8, 6))
        
        self.set_x_ticks_and_labels(ax, pKD_begin, pKD_end)
        
        for idx in reversed(range(inhibitor_kds.shape[0])):
            ax.plot(x_axis, y[idx], 'k', label=self.float_to_prettyprint_conc(kdpi[idx]),
                    marker=self.plot_marker_styles[idx], markevery=self.x_axis_resolution//20, linewidth=1)
        ax.set_xlabel(r"Ligand pK$_\mathrm{D}$")
        ax.set_ylabel("Fraction ligand bound")
        ax.legend()
        ax.grid()
        ax.title.set_text(f"Fraction ligand bound over a range of ligand {self.kd_str}s, [L$_0$]={self.float_to_prettyprint_conc(l)}, [I$_0$]={self.float_to_prettyprint_conc(i)}" +
                          f"\nTarget fraction ligand bound without inhibitor = {target_fraction_ligand_bound}")
        ax.set_xlim(3, 12)
        ax.vlines(6.975,0,1,linestyles="--")
        ax.set_ylim(0, target_fraction_ligand_bound*1.1)
        plt.show()

    def plot_ligand_kd_vs_FLB_as_percentage(
        self,
        pKD_begin: float = 3,
        pKD_end: float = 12,
        l: float = 10e-9,
        i: float = 10e-6,
        kdpi: Union[float, List[float], Tuple[float]] = [1e-9, 10e-9, 100e-9, 1e-6],
        target_fraction_ligand_bound=0.7,
    ):
        if isinstance(kdpi, (float, int)):
            kdpi=[kdpi]
        x_axis = np.linspace(pKD_begin, pKD_end, self.x_axis_resolution)
        ligand_kd_range = 10**(-x_axis)  # We are working in µM, which is 1e-6.
        kdpi = np.array(kdpi)
        protein_concs = np.array(calc_amount_p(
            target_fraction_ligand_bound,l, ligand_kd_range))
        y = np.full((kdpi.shape[0], x_axis.shape[0]), np.nan)


        for it_inhibitor_kds, inhibitor_kd in enumerate(kdpi):
            for it_ligand_kd_range, ligand_kd in enumerate(ligand_kd_range):
                p = protein_concs[it_ligand_kd_range]
                y[it_inhibitor_kds][it_ligand_kd_range] = competition_pl(
                    **{'p': p, 'l': l, 'i': i, 'kdpl': ligand_kd, 'kdpi': inhibitor_kd})/l
        y = ((target_fraction_ligand_bound-y)/target_fraction_ligand_bound)*100

        fig, ax = plt.subplots(figsize=(8, 6))
        self.set_x_ticks_and_labels(ax, pKD_begin, pKD_end)
        plot_line_labels = self.get_plot_line_labels(kdpi)
        for idx in range(kdpi.shape[0]):
            ax.plot(x_axis, y[idx], 'k', label=plot_line_labels[idx],
                    marker=self.plot_marker_styles[idx], markevery=self.x_axis_resolution//10+1, linewidth=1)
        #   ax.plot(x_axis, protein_concs/10, label="[P]")
        ax.set_xlabel(f"Ligand p{self.kd_str}")
        ax.set_ylabel("% Signal")
        ax.legend()
        ax.title.set_text(f"Protein-ligand signal over a range of ligand {self.kd_str}s, [L]={self.float_to_prettyprint_conc(l)}, [I]={self.float_to_prettyprint_conc(i)}"
            +f"\nTarget fraction ligand bound without inhibitor = {target_fraction_ligand_bound}")
        ax.set_xlim(pKD_begin, pKD_end)
        ax.set_ylim(0, 100*1.025)
        plt.show()
    