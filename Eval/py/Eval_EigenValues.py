import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
import csv

FileName = "EigenValues_M14450_Tres255"
Path = "/home/jakob/Downloads/TST_MKL_Eigen/TST/cmake-build-release-intel-2019/"
with open(Path + FileName + '.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    x, y = [], []
    for index_1, row in enumerate(spamreader):
        if index_1 > 1:
            t = float(row[0])
            for element in row[1:]:
                elm = float(element)
                if 10 ** -8 > elm and elm > 0:
                    y.append(elm)
                    x.append(t)
        else:
            print(row)

y = y - np.sum(y) / len(y)
plt.plot(x[:int(len(x)/2)], np.abs(fft(y)[:int(len(y)/2)]))
#plt.plot(x, )
# plt.xscale("log")
# plt.ylim(-1, 1)
# plt.yscale('log')
# plt.xlim(0,1.5)
# ax.set_rmax()
plt.show()
# plt.savefig(FileName + "jpg", quality=100, optimize=True, dpi=600)
