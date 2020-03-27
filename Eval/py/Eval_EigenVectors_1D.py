import csv
from colorsys import hls_to_rgb
from typing import Iterable, Tuple, TypeVar

import numpy as np
import pylab as plt
from numpy import pi

T = TypeVar("T")


def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)


def colorize(z):
    r = np.abs(z)
    arg = np.angle(z)

    h = 180 * arg / pi  # (arg + pi) / (2 * pi) + 0.5
    l = 1.0 - 1.0 / (1.0 + r ** 0.3)
    s = 0.8

    c = np.vectorize(hls_to_rgb)(h, l, s)  # --> tuple
    #c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    #c = c.swapaxes(0, 2)
    return c


FileName = "EigenVectors_M202_Tres400"
Path = "/home/jakob/CLionProjects/TST_MKL_Eigen/TST/cmake-build-release/"
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

b1,b2 = [],[]
l = len(a[0])
l2 = int(len(a[0]) / 2)

print("l=",l)

plt.savefig("KitaevMajoranas+" + ".png", quality=100, optimize=True, dpi=600)
for i in range(int(len(a))):
    print(i)
    b = a[i*2]
    img = []
    for k in range(l2):
        x = b[k]
        y = b[k+l2]
        img.append(x.__abs__()+y.__abs__())

    # 'nearest' interpolation - faithful but blocky
    plt.plot(img) #, interpolation='none', aspect=int(l2/4))
    plt.savefig("Kitaev/Modes" + str(i) + ".png", quality=90, optimize=True)
    plt.clf()
    del b, k, img
