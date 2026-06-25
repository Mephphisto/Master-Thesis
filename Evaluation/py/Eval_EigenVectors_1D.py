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
Path = "/Users/jakobteuffel/Master-Thesis/Simulation/build/"
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
idx = 0
lower, upper = 480,520
for i in range(int(len(a))):
    print(i)
    b = np.array(a[i])
    xs, ys = np.empty(l2, complex), np.empty(l2,complex)
    if (np.sum(np.abs(b[lower: upper])) > .1):
        for k in range(l2):
            x = b[k]
            y = b[k + l2]
            xs[k] = x
            ys[k] = y

        # 'nearest' interpolation - faithful but blocky
        #plt.ylim(0, 0.1)
        #plt.xlim(lower, upper)
        plt.ylim(0,1)
        plt.xlim(lower,upper)
        plt.plot(np.abs(xs+ys))
        plt.plot(np.abs(xs-ys))
        plt.savefig("Kitaev/Modes" + str(idx) + ".png", quality=90, optimize=True)
        plt.clf()
        plt.plot(np.abs(fft(xs-ys))[:l4])
        plt.plot(np.abs(fft(xs+ys)[:l4]))
        plt.ylim(0,13)
        plt.savefig("Kitaev_FFT/Modes_FFT" + str(idx) + ".png", quality=90, optimize=True)
        print("[", idx, "]")
        idx = idx + 1
        plt.clf()
        del k, xs, ys
