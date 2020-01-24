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

    h = 180 * arg / pi  # (arg + pi) / (2 * pi) + 0.5
    l = 1.0 - 1.0 / (1.0 + r ** 0.3)
    s = 0.8
    l = l*4
    c = np.vectorize(hls_to_rgb)(h, l, s)  # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    c = c.swapaxes(0, 2)
    return c


FileName = "EigenVectors_M800_Tres1"
Path = "/Users/jakobteuffel/CodeBase/TST_MKL_Eigen/TST/cmake-build-debug/"
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

l = len(a)
l2 = int(len(a) / 2)
gsize = int(math.sqrt(l2))
l4 = int(len(a) / 4)
print()
for i in range(l):
    b1 = []
    for k in range(l):
        b1.append(a[i][k])
    img = []
    for k in range(gsize):
        aux = []
        for j in range(gsize):
            aux.append(b1[k + gsize * j] * b1[k + gsize * j + l2])
        img.append(aux)
    # 'nearest' interpolation - faithful but blocky
    plt.figure(figsize=(gsize, gsize))
    plt.imshow(colorize(img), interpolation='none')
    plt.savefig("Modes/2D_Modes_" + str(i) + ".png", quality=90, optimize=True)
    plt.clf()
    sys.stdout.write("Progress= " + str(i / l * 100) + "%\r")

    del b1, k, img, aux
