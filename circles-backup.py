import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, RegularPolygon, Ellipse
import numpy as np

# (modified from one of the matplotlib gallery examples)
resolution = 50 # the number of vertices
N = 100
x       = np.random.random(250)
y       = np.random.random(72)
radii   = 0.1*np.random.random(250*72)
patches = []

red = (1.0, 0.0 , 0.0)

print red
print type(red)

red = red + (0.5,)

print red


for x1,y1,r in zip(x, y, radii):
	print (x1,y1)
	circle = Ellipse((x1,y1), r, 0.1)
	print circle
	patches.append(circle)

fig = plt.figure()
ax = plt.axes()
#ax = fig.add_subplot(111)

colors = 100*np.random.random(250*72)
p = PatchCollection(patches, cmap=matplotlib.cm.rainbow, alpha=0.5)
p.set_array(colors)
ax.add_collection(p)
plt.colorbar(p)

plt.show()