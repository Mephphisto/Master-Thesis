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

img = np.empty((100,100),complex)
for x in range(100):
    for y in range(100):
        img[x,y] = (complex((x-50)/5, (y-50)/5))

plt.imshow(colorize(img), interpolation='none',extent=[-10,10,-10,10])
plt.xlabel("Real")
plt.ylabel("Imag")
plt.savefig("Color Wheel", quality=90, optimize=True)
plt.clf()