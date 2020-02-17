import csv
from colorsys import hls_to_rgb
from typing import Iterable, Tuple, TypeVar
import math
import numpy as np
import pylab as plt
from numpy import pi
import sys

T = TypeVar("T")


def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)


def colorize(z):
    r = np.abs(z)
    arg = np.angle(z)

    h = (arg + pi)  / (2 * pi) + 0.5
    l = 1.0 - 1.0/(1.0 + r**0.3)
    s = 0.8

    c = np.vectorize(hls_to_rgb) (h,l,s) # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    c = c.swapaxes(0,2)
    return c


FileName = "DELTA_M1922_Tres32"
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

plt.figure(figsize=(gsize, gsize))

majoranas = range(len(a))

for i in majoranas:
    print(i)
    b1 = []
    for k in range(l):
        b1.append(a[i][k])
    img = []
    for j in range(gsize):
        aux = []
        for k in range(gsize):
            x = b1[k + gsize * j]
            y = b1[k + gsize * j + l2]
            aux.append((x-y))
        img.append(aux)
    # 'nearest' interpolation - faithful but blocky
    plt.imshow(colorize(img), interpolation='none')
    plt.savefig("Majoranas/2D_Modes_" + str(i) + ".png", quality=90, optimize=True)
    plt.clf()
    del b1, k, img, aux

for i in majoranas:
    print(i)
    b1 = []
    for k in range(l):
        b1.append(a[i][k])
    img = []
    for j in range(gsize):
        aux = []
        for k in range(gsize):
            x = b1[k + gsize * j]
            y = b1[k + gsize * j + l2]
            aux.append(y-x)
        img.append(aux)
    # 'nearest' interpolation - faithful but blocky
    plt.imshow(colorize(img), interpolation='none')
    plt.savefig("Majoranas/2D_Modes_-" + str(i) + ".png", quality=90, optimize=True)
    plt.clf()
    del b1, k, img, aux



for i in majoranas:
    print(i)
    b1 = []
    for k in range(l):
        b1.append(a[i][k])
    img = []
    for j in range(gsize):
        aux = []
        for k in range(gsize):
            x = b1[k + gsize * j]
            y = b1[k + gsize * j + l2]
            aux.append(x*x.conjugate()+y.conjugate()*y)
        img.append(aux)
    # 'nearest' interpolation - faithful but blocky
    plt.imshow(colorize(img), interpolation='none')
    plt.savefig("Majoranas/2D_Modes_A" + str(i) + "Density.png", quality=90, optimize=True)
    plt.clf()
    del b1, k, img, aux