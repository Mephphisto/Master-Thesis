import csv
from colorsys import hls_to_rgb
from typing import Iterable, Tuple, TypeVar

import numpy as np
import pylab as plt
from scipy.fft import fft
from numpy import pi

T = TypeVar("T")
def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)


FileName = "EigenVectors_M2002_Tres800"
Path = "/home/jakob/CLionProjects/TST_MKL_Eigen/TST/cmake-build-release-intel/"
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

b1, b2 = [], []
l = len(a[0])
l2 = int(len(a[0]) / 2)
l4 = int(l/4)

print("l=", l)
b = np.array(a[0])
plt.savefig("KitaevMajoranas+" + ".png", quality=100, optimize=True, dpi=600)
lower, upper = 480,520
b_last = np.empty(l)
first = True
xs,ys = [], []
for i in range(int(len(a))):
    print(i)
    b = np.array(a[i])
    b = b / np.linalg.norm(b)
    if (np.sum(np.abs(b[lower: upper])) > .1):
        b_last = b
        if first:
            first = False
            ys.append(1)
        else:
            x = np.vdot(b, b_last).__abs__()
            print(x)
            xs.append(x)
            ys.append(ys[-1]*x)

plt.plot(xs)
plt.plot(ys)
plt.yscale("log")
print(ys[-1])
plt.show()
