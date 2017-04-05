import sys
 
import numpy as np
from numpy import array as npa
 
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, CircleCollection
from matplotlib.gridspec import GridSpec

from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, RegularPolygon, Ellipse, CirclePolygon
 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital


bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
data = Procar("./PROCAR").data
font = {'family': 'serif', 'size': 48}
plt.rc('font', **font)

data = data[Spin.up].sum(axis=2)


N = 72*250
k       = np.arange(250)

radii   = 0.1*np.random.random(250*72)

patches = []


for b in range (72):
	listek = zip(range(len(bands.kpoints)), bands.bands[Spin.up][b] - bands.efermi)
	
	#print listek[1][1]
	for k in range(250):
		circle = Circle(listek[k], 0.01)
		print circle
		patches.append(circle)

fig = plt.figure()
ax = plt.axes()

ax.set_xlim(1, 250)
ax.set_ylim(-4.0, 4.0)

colors = 100*np.random.random(250*72)
p = PatchCollection(patches, cmap=mpl.cm.rainbow, alpha=0.8)
p.set_array(colors)
ax.add_collection(p)
plt.colorbar(p)
plt.show()