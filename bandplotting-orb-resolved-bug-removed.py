#!/usr/bin/env python
# -*- coding=utf-8 -*-
#   A Python code for plotting orbital-resolved bandstructure. 
#   Written by              : Internet
#   catalyst                : Ravindra
#   under the eagle eyes of : Rinkle

import sys
import os 

import numpy as np
from numpy import array as npa
 
import matplotlib as mpl
import matplotlib.pyplot as plt

 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital


mpl.rc('text', usetex=True)
mpl.rc('font', weight='bold')
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
 
 
if __name__ == "__main__":
    # Data is read from vasprun.xml and PROCAR files using pymatgen library
    # ---------
 
    # kpoints labels can be read from the built-in functions. Right now, set manually.
    # path = HighSymmKpath(mg.Structure.from_file("./CONTCAR")).kpath["path"]

    path = [["X", "\Gamma", "R", "A", "Z", "\Gamma"]]
    labels = ["X", r"$\boldsymbol{\Gamma}$", "R", "A", "Z", r"$\boldsymbol{\Gamma}$"]  #[r"$%s$" % lab for lab in path[0][0:6]]
 
    # bands object prepared using pymatgen library. contains eigenvalue information
    bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
 
    # projected bands read from procar files
    data = Procar("./PROCAR").data
 
    # set up matplotlib plot
    # ----------------------
 
    # general options for plot
    font = {'family': 'serif', 'size': 20, 'weight': 'bold'}
    plt.rc('font', **font)
     
    markers     = ("o", "^", "h", "*", "s", "D")
    colorlist   = ("r", "g", "b", "m", "c", "y")

    # set ylim for the plot
    # ---------------------
    ax = plt.axes()
 
    emin = -1 
    emax =  1 
    ax.set_ylim(emin, emax)

    # Band Diagram
    # ------------
 
    # sum up contribution over carbon atoms
    data = data[Spin.up].sum(axis=2)
 
    # sum up px and py contributions and normalize contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 4))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = data[k][b][Orbital.s.value]**2
            pxc = data[k][b][Orbital.px.value]**2 
            pyc = data[k][b][Orbital.py.value]**2
            pzc = data[k][b][Orbital.pz.value]**2
            tot = sc + pxc + pyc + pzc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pxc / tot
                contrib[b, k, 2] = pyc / tot
                contrib[b, k, 3] = pzc / tot
                #print k, b, sum(contrib[b,k,:])
                

    #colors  = 0.5*np.ones(len(bands.kpoints)*bands.nb_bands)
    patches = []
    

    # for b in range(bands.nb_bands):
    #     coords = zip(range(len(bands.kpoints)), bands.bands[Spin.up][b] - bands.efermi)
    #     for k in range(len(bands.kpoints)):
    #         x = coords[k][0]
    #         y = coords[k][1]
    #         index_max = np.argmax(np.array([contrib[b, k, 0], contrib[b, k, 1], contrib[b, k, 2], contrib[b, k, 3]]))
    #         plt.scatter(x,y, s = 50*contrib[b,k,index_max], marker = markers[index_max], color=colorlist[index_max], alpha=0.4)

    #fig = plt.figure()



    for b in range(bands.nb_bands):
        coords = zip(range(len(bands.kpoints)), bands.bands[Spin.up][b] - bands.efermi)
        for k in range(len(bands.kpoints)):
            x = coords[k][0]
            y = coords[k][1]
            #index_max = np.argmax(np.array([contrib[b, k, 0], contrib[b, k, 1], contrib[b, k, 2], contrib[b, k, 3]]))
            plt.scatter(x,y, s = 50*contrib[b,k,0], marker = markers[0], color=colorlist[0], alpha=0.3)
            plt.scatter(x,y, s = 50*contrib[b,k,1], marker = markers[1], color=colorlist[1], alpha=0.3)
            plt.scatter(x,y, s = 50*contrib[b,k,2], marker = markers[2], color=colorlist[2], alpha=0.3)
            plt.scatter(x,y, s = 50*contrib[b,k,3], marker = markers[3], color=colorlist[3], alpha=0.3)



    ax.set_ylabel(r"E - E$_f$ (eV)", fontsize=20 , weight='bold')
    ax.set_xticklabels(['-1.0', '-0.5', '0.0', '0.5', '1.0'], fontsize=20, weight='bold')
    ax.yaxis.label.set_fontsize(20) 

    ax.grid()
 
    # fermi level at 0
    ax.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="0.75", lw=0.5)
 
    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        plt.vlines(i * step, emin, emax, "0.75")
    ax.set_xticks([i * step for i in range(nlabs)])
    ax.set_xticklabels(labels)
    ax.xaxis.label.set_fontsize(20) 
 
    ax.set_xlim(0, len(bands.kpoints))
 
    ## Following part is for hiding the x and y ticklabels
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)

    #graph_title = os.path.basename(os.path.dirname(os.path.realpath('vasprun.xml'))) 
    #direct = str(os.getcwd())
    #graph_title = os.system(direct.replace("/", "-"))
    #ax.set_title("Ge-based--As2GeCd--band-PBE")
    #plt.show()
    plt.savefig(sys.argv[1] + ".pdf", format="pdf")
