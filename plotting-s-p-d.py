#!/usr/bin/env python
# -*- coding=utf-8 -*-
 
import sys
 
import numpy as np
from numpy import array as npa
 
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
 
 
def rgbline(ax, k, e, red, green, blue, alpha=1.):
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
 
    nseg = len(k) - 1

    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linestyle = "solid", linewidth=8)
    ax.add_collection(lc)
 
if __name__ == "__main__":
    # kpoints labels
    # path = HighSymmKpath(mg.Structure.from_file("./CONTCAR")).kpath["path"]
    path = [["X", "\Gamma", "R", "A", "Z", "\Gamma"]]
    labels = [r"$%s$" % lab for lab in path[0][0:6]]
 
    # bands
    bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
 
    # projected bands
    data = Procar("./PROCAR").data
 
    # set up matplotlib plot
    # ----------------------
 
    # general options for plot
    font = {'family': 'serif', 'size': 48}
    plt.rc('font', **font)
 
    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of State
    gs = GridSpec(1, 1)
    fig = plt.figure(figsize=(23, 15))
    fig.suptitle("Sb$_2$SnZn (no SOC)")
    ax1 = plt.subplot(gs[0])
    #ax2 = plt.subplot(gs[1])  # , sharey=ax1)
 
    # set ylim for the plot
    # ---------------------
    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))
 
    emin -= bands.efermi + 1
    emax -= bands.efermi - 1
    ax1.set_ylim(-1, 1)
    #ax2.set_ylim(emin, emax)
 
    # Band Diagram
    # ------------
 
    # sum up contribution over carbon atoms
    data = data[Spin.up].sum(axis=2)
 
    # sum up px and py contributions and normalize contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = data[k][b][Orbital.s.value]**2
            pc = data[k][b][Orbital.px.value]**2 + data[k][b][Orbital.py.value]**2 + data[k][b][Orbital.pz.value]**2 
            dc = data[k][b][Orbital.dxy.value]**2 + data[k][b][Orbital.dxy.value]**2 + data[k][b][Orbital.dyz.value]**2 + data[k][b][Orbital.dxz.value]**2 + data[k][b][Orbital.dz2.value]**2 + data[k][b][Orbital.dx2.value]**2 

            tot = sc + pc + dc 
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pc / tot
                contrib[b, k, 2] = dc / tot
 
 
    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1,
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0],
                contrib[b, :, 1],
                contrib[b, :, 2],)
 
    # style

    ax1.set_ylabel(r"$E - E_f$ (eV)")
    ax1.grid()
    

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="0.75", lw=0.5)
    #ax2.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=0.5)
 
    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        ax1.vlines(i * step, emin, emax, color= "0.75", lw=0.5)
    ax1.set_xticks([i * step for i in range(nlabs)])
    ax1.set_xticklabels(labels)
 
    ax1.set_xlim(0, len(bands.kpoints))
 
    ax1.text(0.9, 0.62, '__ $p_x$', verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='red', fontsize=48)
    ax1.text(0.9, 0.57, '__ $p_y$', verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='green', fontsize=48)
    ax1.text(0.9, 0.52, '__ $p_z$', verticalalignment='bottom', horizontalalignment='left', transform=ax1.transAxes, color='blue', fontsize=48)
    # plt.show()
    
    plt.savefig(sys.argv[1] + ".eps", format="eps")
    

