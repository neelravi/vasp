#!/usr/bin/env python
# -*- coding=utf-8 -*-
 
import sys
 
import numpy as np
from numpy import array as npa
 
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, CircleCollection
from matplotlib.gridspec import GridSpec

from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, RegularPolygon, Ellipse
 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
 
 
def rgbline1(ax, k, e, red, green, blue, alpha=0.5):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
 
    nseg = len(k) - 1

    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=8)
    ax.add_collection(lc)

#def patchup(ax, k, e, red, green, blue, alpha=0.5):
def patchup(ax, k, e):
    patches = []

    radii   = 0.1*np.random.random(250)

    for i in range(250):
        circle = Circle((i,e), radii[i])
        patches.append(circle)

    colors = 250*np.random.random(250)
    #colorred = [100 * (red[i]) for i in range(250)]
    p = PatchCollection(patches,  cmap=mpl.cm.jet , alpha=0.4)
    p.set_array(colors)
    ax.add_collection(p)
    plt.colorbar(p)
 #   g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
 #   b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
 #   a = np.ones(nseg, np.float) * alpha
    #cc = CircleCollection(redarea)
    #cc.set_color('r')
    #ax.add_collection(cc, autolim=True)


 
if __name__ == "__main__":
    # read data
    # ---------
 
    # kpoints labels
#    path = HighSymmKpath(mg.Structure.from_file("./CONTCAR")).kpath["path"]
    path = [["X", "\Gamma", "R", "A", "Z", "\Gamma"]]
    labels = [r"$%s$" % lab for lab in path[0][0:6]]
 
    # bands
    bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
 
    # projected bands
    data = Procar("./PROCAR").data
 
    # density of state
   # dosrun = Vasprun("./vasprun.xml")
 
    # set up matplotlib plot
    # ----------------------
 
    # general options for plot
    font = {'family': 'serif', 'size': 48}
    plt.rc('font', **font)
 
    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of State
    #gs = GridSpec(1, 1, width_ratios=[1, 1])
    #fig = plt.figure(figsize=(23, 15))
#    fig.suptitle("Bands diagram of As2SnZn")
    #ax1 = plt.subplot(gs[0])
#    ax2 = plt.subplot(gs[1])  # , sharey=ax1)
 
    # set ylim for the plot
    # ---------------------
    ax = plt.axes()
    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))
 
    emin -= bands.efermi + 1
    emax -= bands.efermi - 1
    ax.set_ylim(-4, 4)
    #ax2.set_ylim(emin, emax)
 
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
 
    # plot bands using rgb mapping

    radii   = 0.1*np.random.random(250*72)
    colors  = 250*72*np.random.random(250*72)
    patches = []
    print len(radii), len(colors), len(patches)
    print bands.efermi
    tup1 = None

    for b in range(bands.nb_bands):
        coords = zip(range(len(bands.kpoints)), [e - bands.efermi for e in bands.bands[Spin.up][b]])

        for k in range (250):
            circle = Ellipse(coords[k], 2.0, 0.1)
            #print circle
            patches.append(circle)  
    

    fig = plt.figure()


    p = PatchCollection(patches,mpl.cm.rainbow)     
    p.set_array(colors)
    ax.add_collection(p)
    fig.colorbar(p,ax=ax)
#    print zip(kpts , energ)
    

#    p = PatchCollection(patches, cmap=mpl.cm.rainbow)
#    p.set_array(colors)
#    ax.add_collection(p)
#    plt.colorbar(p)



    # for b in range(bands.nb_bands):

    #     for kpts, energ, rad in  zip(range(len(bands.kpoints)), [e - bands.efermi for e in bands.bands[Spin.up][b]], radii):

    #         circle = Circle((kpts, energ), radii)
    #         patches.append(circle)
            
    #     p = PatchCollection(patches,  cmap=mpl.cm.jet , alpha=0.4)
    #     p.set_array(colors)
    #     ax.add_collection(p)
    #     plt.colorbar(p)
            
        #colorred = [100 * (red[i]) for i in range(250)]

        #patchup(ax, range(len(bands.kpoints)), [e - bands.efermi for e in bands.bands[Spin.up][b]],
        #        contrib[b, :, 0], contrib[b, :, 1], contrib[b, :, 2])
 
            #patchup(ax, range(len(bands.kpoints)), [e - bands.efermi for e in bands.bands[Spin.up][b]])

    # style
#    plt.set_xlabel("Ravindra's own kpoints ")
    ax.set_ylabel(r"$E - E_f$ (eV)")
    ax.grid()
 
    # fermi level at 0
    ax.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=0.5)
 
    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        plt.vlines(i * step, emin, emax, "k")
    ax.set_xticks([i * step for i in range(nlabs)])
    ax.set_xticklabels(labels)
 
    ax.set_xlim(0, len(bands.kpoints))
 
    plt.show()
