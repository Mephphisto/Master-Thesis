import csv
from colorsys import hls_to_rgb
from typing import Iterable, Tuple, TypeVar
import math
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
import sys
from mpl_toolkits import mplot3d

T = TypeVar("T")

def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)


FileName = "MU_M1922_Tres32"
Path = "/home/jakob/Downloads/TST_MKL_Eigen/TST/cmake-build-release-intel-2019/"
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


X, Y = np.meshgrid(np.linspace(0, gsize - 1, gsize), np.linspace(0, gsize - 1, gsize))

for i in range(len(a)):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    print(i)
    b1 = []
    for k in range(l):
        b1.append(a[i][k])
    Z = np.empty((gsize, gsize))
    for j in range(gsize):
        for k in range(gsize):
            x = b1[k + gsize * j]
            y = b1[k + gsize * j + l2]
            z = (x.real - y.real)/2
            Z[k, j] = float(z.real)

    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Density')
    ax.view_init(20, 25)
    fig.show()
    fig.savefig("Majoranas2D/2D_Modes_" + str(i) + ".png", quality=90, optimize=True)
    fig.clf()
    ax.clear()
    del b1, Z
