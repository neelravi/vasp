#import sys
 
import numpy as np
from numpy import array as npa
 
import matplotlib.pyplot as plt
#from matplotlib.collections import LineCollection
#from matplotlib.gridspec import GridSpec
 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital


#bands = Vasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
dosrun = Vasprun("./vasprun.xml")
spd_dos = dosrun.complete_dos.get_spd_dos()


run = Vasprun("./vasprun.xml", parse_projected_eigen=True)
bands = run.get_band_structure("./KPOINTS", line_mode=True, efermi=dosrun.efermi)


data = Procar("./PROCAR").data
data = data[Spin.up].sum(axis=2)


klist = [[i] for i in range(250)]

name = [None]*3


name[0] = "As"
name[1] = "Ge"
name[2] = "Cd"

##pbands = bands.get_projections_on_elts_and_orbitals({name: ["s", "p", "d"]})
pbands0 = bands.get_projections_on_elements_and_orbitals({name[0]: ["s", "p"]})
pbands1 = bands.get_projections_on_elements_and_orbitals({name[1]: ["s", "p"]})
pbands2 = bands.get_projections_on_elements_and_orbitals({name[2]: ["s", "p"]})

# compute s, p, d normalized contributions
contrib = np.zeros((3, bands.nb_bands, len(bands.kpoints), 2))
for b in range(bands.nb_bands):
    for k in range(len(bands.kpoints)):
        	sc0 = pbands0[Spin.up][b][k][name[0]]["s"]**2
        	pc0 = pbands0[Spin.up][b][k][name[0]]["p"]**2
        	tot0 = sc0 + pc0
        	if tot0 != 0.0:
				contrib[0, b, k, 0] = sc0 / tot0
				contrib[0, b, k, 1] = pc0 / tot0

        	sc1 = pbands1[Spin.up][b][k][name[1]]["s"]**2
        	pc1 = pbands1[Spin.up][b][k][name[1]]["p"]**2
        	tot1 = sc1 + pc1
        	if tot1 != 0.0:
				contrib[1, b, k, 0] = sc1 / tot1
				contrib[1, b, k, 1] = pc1 / tot1

        	sc2 = pbands2[Spin.up][b][k][name[2]]["s"]**2
        	pc2 = pbands2[Spin.up][b][k][name[2]]["p"]**2
        	tot2 = sc2 + pc2
        	if tot2 != 0.0:
				contrib[2, b, k, 0] = sc2 / tot2
				contrib[2, b, k, 1] = pc2 / tot2


markers     = ("o", "^", "s", "h", "D", "*")
colorlist   = ("r", "g", "b", "m", "c", "y")

for j in range(72):
	ek = zip(range(len(bands.kpoints)), [e - bands.efermi for e in bands.bands[Spin.up][j]]) 

	orbital_contribution = zip ([contrib[0,j,k,0] for k in range (250)], [contrib[0,j,k,1] for k in range (250)])
	
	if contrib[0,j,k,0] > contrib[0,j,k,1]:
		marker = "o" ; c = "y"   # s orbital dominating
	else:
		marker = "^" ; c = "m"	# p orbital dominating
	plt.scatter(*zip(*ek),  s = 100, marker = marker, c=c, alpha=0.5)

	# for k in range(len(bands.kpoints)):
	# 	if contrib[1,j,k,0] > contrib[1,j,k,1]:
	# 		marker = "s" ; c = "y"	# s orbital dominating
	# 	else:
	# 		marker = "*" ; c = "m"
	# 	plt.scatter((k,bands.bands[Spin.up][j]-bands.efermi),  s = 100, marker = marker, c=c, alpha=0.5)

	# for k in range(len(bands.kpoints)):
	# 	if contrib[2,j,k,0] > contrib[2,j,k,1]:
	# 		marker = "s" ; c = "y"	# s orbital dominating
	# 	else:
	# 		marker = "*" ; c = "m"
	# 	plt.scatter((k,bands.bands[Spin.up][j]-bands.efermi),  s = 100, marker = marker, c=c, alpha=0.5)
#		plt.scatter(*zip(*ek),  s = 100, marker = marker, c=c, alpha=0.5)

plt.show()