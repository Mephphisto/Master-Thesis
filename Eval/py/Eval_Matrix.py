import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from typing import Iterable, Tuple, TypeVar
from numpy import pi
import pylab as plt
from colorsys import hls_to_rgb

T = TypeVar("T")


def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)


def colorize(z):
    r = np.abs(z)
    arg = np.angle(z)

    h = (arg + pi) / (2 * pi) + 0.5
    l = 1.0 - 1.0 / (1.0 + r ** 0.3)
    s = 0.8

    c = np.vectorize(hls_to_rgb)(h, l, s)  # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    c = c.swapaxes(0, 2)
    return c


FileName = "Matrix320_Tres64"
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

# 'nearest' interpolation - faithful but blocky
plt.imshow(colorize(a), interpolation='nearest')

plt.show(dpi=1200)
# plt.savefig(FileName + ".jpg", quality=100, optimize=True, dpi=600)
