#!/usr/bin/env python
# -*- coding=utf-8 -*-
 
import sys
 
import numpy as np
from numpy import array as npa

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
 

mpl.rc('text', usetex=True)
mpl.rc('font', weight='bold')
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


def rgbline(ax, k, e, red, green, blue, alpha=1.):
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
 
    nseg = len(k) - 1

    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg , linewidth = 4 , colors=list(zip(r, g, b, a)))
    ax.add_collection(lc)


 
if __name__ == "__main__":
    # read data
    # ---------
 
    # kpoints labels
#    path = HighSymmKpath(mg.Structure.from_file("./CONTCAR")).kpath["path"]
    path = [["X", "\Gamma", "R", "A", "Z", "\Gamma"]]
    labels = ["X", r"$\boldsymbol{\Gamma}$", "R", "A", "Z", r"$\boldsymbol{\Gamma}$"]  #[r"$%s$" % lab for lab in path[0][0:6]]
 
    # bands
    bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
 
    # projected bands
    data = Procar("./PROCAR").data
 
    # density of state
    # dosrun = Vasprun("./vasprun.xml")
 
    # set up matplotlib plot
    # ----------------------
 
    # general options for plot
    font = {'family': 'serif', 'size': 20}
    plt.rc('font', **font)
 
    ax = plt.axes()
    # set ylim for the plot
    # ---------------------
    emin = -1 
    emax =  1 
    ax.set_ylim(emin, emax)
 
    # Band Diagram
    # ------------
 
    # sum up contribution over carbon atoms
    data = data[Spin.up].sum(axis=2)
 
    # sum up px and py contributions and normalize contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            pxc = data[k][b][Orbital.px.value]**2 
            pyc = data[k][b][Orbital.py.value]**2
            pzc = data[k][b][Orbital.pz.value]**2
            tot = pxc + pyc + pzc
            if tot != 0.0:
                contrib[b, k, 0] = pxc / tot
                contrib[b, k, 1] = pyc / tot
                contrib[b, k, 2] = pzc / tot
 
    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax,
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0],
                contrib[b, :, 1],
        		contrib[b, :, 2])
 
    # style
#    ax1.set_xlabel("Ravindra's own kpoints ")
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
        ax.vlines(i * step, emin, emax, "0.75", lw=0.25)
    ax.set_xticks([i * step for i in range(nlabs)])
    ax.set_xticklabels(labels)
 
    ax.set_xlim(0, len(bands.kpoints))

    plt.savefig(sys.argv[1] + ".pdf", format="pdf")

