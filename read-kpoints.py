#!/usr/bin/env python
# -*- coding=utf-8 -*-
 
import sys
 
import numpy as np
from numpy import array as npa
 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
 


 
if __name__ == "__main__":
    # read data
    # ---------
 
    # kpoints labels
#    path = HighSymmKpath(mg.Structure.from_file("./CONTCAR")).kpath["path"]
#    path = [["X", "\Gamma", "R", "A", "Z", "\Gamma"]]
#    labels = ["X", r"$\boldsymbol{\Gamma}$", "R", "A", "Z", r"$\boldsymbol{\Gamma}$"]  #[r"$%s$" % lab for lab in path[0][0:6]]
 
    # bands
    #bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)

    a = Vasprun("./vasprun.xml").eigenvalue_band_properties()
    #print a[:]
 
    # projected bands
    # data = Procar("./PROCAR").data
 
    # density of state
    # dosrun = Vasprun("./vasprun.xml")
 
 
    # sum up contribution over carbon atoms
    # data = data[Spin.up].sum(axis=2)
 
    # # sum up px and py contributions and normalize contributions
    # contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    # for b in range(bands.nb_bands):
    #     for k in range(len(bands.kpoints)):
    #         pxc = data[k][b][Orbital.px.value]**2 
    #         pyc = data[k][b][Orbital.py.value]**2
    #         pzc = data[k][b][Orbital.pz.value]**2
    #         tot = pxc + pyc + pzc
    #         if tot != 0.0:
    #             contrib[b, k, 0] = pxc / tot
    #             contrib[b, k, 1] = pyc / tot
    #             contrib[b, k, 2] = pzc / tot
 
    # # plot bands using rgb mapping
    # for b in range(bands.nb_bands):
    #     rgbline(ax,
    #             range(len(bands.kpoints)),
    #             [e - bands.efermi for e in bands.bands[Spin.up][b]],
    #             contrib[b, :, 0],
    #             contrib[b, :, 1],
    #     		contrib[b, :, 2])
 

