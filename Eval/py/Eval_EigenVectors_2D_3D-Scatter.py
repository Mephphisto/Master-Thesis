import csv
from colorsys import hls_to_rgb
from typing import Iterable, Tuple, TypeVar
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate

T = TypeVar("T")


def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)


FileName = "MU_M5000_Tres1"
Path = "/home/jakob/CLionProjects/TST_MKL_Eigen/TST/cmake-build-release-gcc/"
a = []
with open(Path + FileName + '.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for index_1, row in enumerate(spamreader):
        if index_1 > 1:
            row_real = []
            for r, i in grouped(row[1:]):
                r = r.replace("(", "")
                i = i.replace(")", "")
                row_real.append(complex(float(r), float(i)))
            a.append(row_real)
        else:
            print(row)

l = len(a[0])
l2 = int(len(a[0]) / 2)
gsize = int(math.sqrt(l2))
l4 = int(len(a[0]) / 4)
OvSmpl = 10
X, Y = np.meshgrid(np.linspace(0, gsize - 1, gsize), np.linspace(0, gsize - 1, gsize))

for i in range(len(a)):
    fig = plt.figure()
    ax = plt.axes(projection='3d', xticks=[], yticks=[], zticks=[-1, 0, 1])
    print(i)
    b1 = []
    for k in range(l):
        b1.append(a[i][k])
    Z = np.empty((gsize, gsize))
    for j in range(gsize):
        for k in range(gsize):
            x = b1[k + gsize * j]
            y = b1[k + gsize * j + l2]
            z = x  # (x*x.conjugate() + y*y.conjugate())/2
            Z[k, j] = float(z.real)
        # print(Z.min())
    # xi , yi =  np.mgrid[0:gsize-1:1000j, 0:gsize-1:1000j]
    # zi = interpolate.griddata((X.flatten(), Y.flatten()), Z.flatten(), (xi, yi), method='cubic')
    # surf = ax.plot_surface(xi, yi, zi, cmap='RdBu', edgecolor='None',ccount=1000, rcount=1000, linewidth=0,
    #                       antialiased=False)  # ,cmap='gist_heat' cmap='twilight_shifted'
    surf = ax.plot_surface(X, Y, Z, cmap='RdBu', edgecolor='None', ccount=1000, rcount=1000, linewidth=0,
                           antialiased=False)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('µ')
    # ax.set_zlim(-1, 1)
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    ax.view_init(25, 45)

    fig.show()
    # fig.savefig("Mu_Plot" + str(i) + ".png", dpi=800)
    fig.clf()
    ax.clear()
    del b1, Z
