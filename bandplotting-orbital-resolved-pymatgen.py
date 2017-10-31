#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import os 

import numpy as np
from numpy import array as npa
 
import matplotlib as mpl
import matplotlib.pyplot as plt

 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar, BSVasprun
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.plotter import BSPlotterProjected




mpl.rc('text', usetex=True)
mpl.rc('font', weight='bold')
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
 
 
if __name__ == "__main__":

 
    # bands object prepared using pymatgen library. contains eigenvalue information
    v = BSVasprun("./vasprun.xml",parse_projected_eigen=True)
    bs = v.get_band_structure(line_mode=True)
    #print (bs.is_metal())
    #print (bs.get_band_gap())
    #print (bs.get_direct_band_gap())

    #print (bs.get_projections_on_elements_and_orbitals({'As':['s','p','d']}))
    #promenade = HighSymmKpath.get_kpoints
    #print promenade
    #get_kpoints(bs,line_density=20, coords_are_cartesian=True)
    BSPlotter(bs).show()
    
    #BSPlotter(bs).plot_brillouin()
    #BSPlotter(bs).save_plot(filename="normal-bandstructure.pdf",img_format="pdf",zero_to_efermi=True)

    bsproj = BSPlotterProjected(bs).get_projected_plots_dots_patom_pmorb(dictio={'As':['px','py','pz']}, dictpa={'As':[5,6,7,8]}, sum_atoms={'As':[5,6,7,8]}, sum_morbs={'As':['px','py','pz']})
    bsproj.show()


    # trying new things here

    #bandstruct = BSDOSPlotter(bs_projection="As", dos_projection=None, vb_energy_range=2, cb_energy_range=2, egrid_interval=0.5, rgb_legend=True)

    #BSDOSPlotterProjected(bandstruct)


    # #plt.show()
    # plt.savefig(sys.argv[1] + ".pdf", format="pdf")
