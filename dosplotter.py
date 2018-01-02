#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import os 

import numpy as np
from numpy import array as npa
 
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar, BSVasprun
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import Dos, CompleteDos
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter, BSDOSPlotter
from pymatgen.electronic_structure.plotter import BSPlotterProjected




mpl.rc('text', usetex=True)
mpl.rc('font', weight='bold')
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
 
 
if __name__ == "__main__":

 
    # bands object prepared using pymatgen library. contains eigenvalue information
    dosrun = Vasprun("./vasprun.xml",parse_projected_eigen=True, parse_dos=True)
    bands = Vasprun("./vasprun.xml",parse_projected_eigen=True, parse_dos=True).get_band_structure("./KPOINTS", line_mode=False)
    #bs = dosrun.get_band_structure(line_mode=True)
    #dos = d.get_element_spd_dos()
    #print (bs.is_metal())
    #print (bs.get_band_gap())
    #print (bs.get_direct_band_gap())

    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure()
    ax = plt.subplot(gs[0])
    ax.set_yticklabels([])
    ax.grid()


    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))
 
    emin -= bands.efermi + 1
    emax -= bands.efermi - 1
    ax.set_ylim(emin, emax)
    ax.set_ylim(-0.5, 0.5)


    ax.set_xticks(np.arange(0, 0.8, 0.2))
    ax.set_xticklabels(np.arange(0, 0.8, 0.2))
    ax.set_xlim(1e-6, 0.8)
    ax.hlines(y=0, xmin=0, xmax=1, color="k", lw=2)
    ax.set_xlabel("Density of State")

    #se_contrib_array = np.zeros((2000,1), dtype=np.float)
    #print (np.shape(se_contrib_array))
    #se_contrib_array = npa(dosrun.pdos[Site].sum(axis=2)) 
    #print (se_contrib_array)

    # Contribution of px of 1st type of atom (check the position in POSCAR)
    ax.plot(npa(dosrun.pdos[0][Orbital.px][Spin.up]) + npa(dosrun.pdos[1][Orbital.px][Spin.up]) +  
        npa(dosrun.pdos[2][Orbital.px][Spin.up]) + npa(dosrun.pdos[3][Orbital.px][Spin.up]) + 
        npa(dosrun.pdos[4][Orbital.px][Spin.up]) + npa(dosrun.pdos[5][Orbital.px][Spin.up]) + 
        npa(dosrun.pdos[6][Orbital.px][Spin.up]) + npa(dosrun.pdos[7][Orbital.px][Spin.up]) + 
        npa(dosrun.pdos[8][Orbital.px][Spin.up]) + npa(dosrun.pdos[9][Orbital.px][Spin.up]) + 
        npa(dosrun.pdos[10][Orbital.px][Spin.up]) + npa(dosrun.pdos[11][Orbital.px][Spin.up]),
             dosrun.tdos.energies - dosrun.efermi, "r-", label="Se-px", linewidth=2)
 
    # Contribution of px of 2nd type of atom (check the position in POSCAR)
    ax.plot(npa(dosrun.pdos[12][Orbital.px][Spin.up]) + npa(dosrun.pdos[13][Orbital.px][Spin.up]) +
        npa(dosrun.pdos[14][Orbital.px][Spin.up]) + npa(dosrun.pdos[15][Orbital.px][Spin.up]) +
        npa(dosrun.pdos[16][Orbital.px][Spin.up]) + npa(dosrun.pdos[17][Orbital.px][Spin.up]) +
        npa(dosrun.pdos[18][Orbital.px][Spin.up]) + npa(dosrun.pdos[19][Orbital.px][Spin.up]) +
        npa(dosrun.pdos[20][Orbital.px][Spin.up]) + npa(dosrun.pdos[21][Orbital.px][Spin.up]) +
        npa(dosrun.pdos[22][Orbital.px][Spin.up]) + npa(dosrun.pdos[23][Orbital.px][Spin.up]),
             dosrun.tdos.energies - dosrun.efermi, "g-", label="Te-px", linewidth=2)
 

    # Contribution of px of 3rd type of atom (check the position in POSCAR)
    ax.plot(npa(dosrun.pdos[24][Orbital.px][Spin.up]) + npa(dosrun.pdos[25][Orbital.px][Spin.up]) +  
            npa(dosrun.pdos[26][Orbital.px][Spin.up]) + npa(dosrun.pdos[27][Orbital.px][Spin.up]) +
            npa(dosrun.pdos[28][Orbital.px][Spin.up]) + npa(dosrun.pdos[29][Orbital.px][Spin.up]) +  
            npa(dosrun.pdos[30][Orbital.px][Spin.up]) + npa(dosrun.pdos[31][Orbital.px][Spin.up]) +  
            npa(dosrun.pdos[32][Orbital.px][Spin.up]) + npa(dosrun.pdos[33][Orbital.px][Spin.up]) +  
            npa(dosrun.pdos[34][Orbital.px][Spin.up]) + npa(dosrun.pdos[35][Orbital.px][Spin.up]) , 
             dosrun.tdos.energies - dosrun.efermi, "b-", label="Bi-px", linewidth=2)
    #plt.show()
    plt.savefig(sys.argv[1] + ".pdf", format="pdf") 

    #print (v.pdos[0][Orbital.s][Spin.up], v.pdos[3][Orbital.s][Spin.up], v.pdos[5][Orbital.s][Spin.up])

    #print (bs.get_projections_on_elements_and_orbitals({'As':['s','p','d']}))
    #promenade = HighSymmKpath.get_kpoints
    #print promenade 
    #get_kpoints(bs,line_density=20, coords_are_cartesian=True)
    #BSPlotter(bs).show()
    
    #BSPlotter(bs).plot_brillouin()
    #BSPlotter(bs).save_plot(filename="normal-bandstructure.pdf",img_format="pdf",zero_to_efermi=True)

    #bsproj = BSPlotterProjected(bs).get_projected_plots_dots_patom_pmorb(dictio={'Ge':['s']}, dictpa={'Ge':[4]})
    #bsproj.show()


    # trying new things here

    #bandstruct = BSDOSPlotter(bs_projection="As", dos_projection="As", vb_energy_range=2, cb_energy_range=2, egrid_interval=0.5, rgb_legend=True)
    #print (type(bandstruct))
    #print (bandstruct.get_plot(bs).show())
