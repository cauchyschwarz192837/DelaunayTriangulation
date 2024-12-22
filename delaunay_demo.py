import matplotlib.pyplot as plt
from delaunay import *
from primitives import *

# data source: https://raw.githubusercontent.com/plotly/datasets/master/api_docs/mt_bruno_elevation.csv

data = []

with open('mt_bruno_elevation.dat', 'r') as file:
    for line in file:
        nums = [float(num.strip()) for num in line.split(',')]
        data.append(nums)

xs = []
ys = []
zs = []
pts = []
idx = {}

scale = 10000 # perturbation hack to ensure various assumptions: no 3 collinear, no 2 with same x-coord
for i,row in enumerate(data):
    for j,num in enumerate(row):
        zs.append(num)
        pt = Point(scale*i+j**2,scale*j+i**2,scale)
        xs.append(pt.x())
        ys.append(pt.y())
        pts.append(pt)
        idx[pts[-1]] = len(zs)-1

T = Triangulation(pts, use_tree=False, make_legal=False)
T.random_incremental()
triangles = [ [idx[q] for q in tri.to_tuple()] for tri in T.get_triangles()]

T.draw()
plt.show()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_trisurf(xs, ys, zs, triangles=triangles, cmap=plt.cm.viridis, linewidth=0.2)
plt.show()

# T = Triangulation(pts, use_tree=True, make_legal=True)
# T.random_incremental()
# triangles = [ [idx[q] for q in tri.to_tuple()] for tri in T.get_triangles()]

# T.draw()
# plt.show()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_trisurf(xs, ys, zs, cmap=plt.cm.viridis, linewidth=0.2)
plt.show()