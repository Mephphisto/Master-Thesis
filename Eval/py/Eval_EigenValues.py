import numpy as np
import matplotlib.pyplot as plt
import csv

FileName = "EigenValues_M2002_Tres800"
Path = "/home/jakob/CLionProjects/TST_MKL_Eigen/TST/cmake-build-release-intel/"
with open(Path + FileName + '.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for index_1, row in enumerate(spamreader):
        if index_1 > 1:
            t = float(row[0]),

            x, y = [], []
            for element in row[1:]:
                x.append(t)
                y.append(float(element))

            plt.plot(x, y, marker=".", linewidth=0, markersize=.1)

        else:
            print(row)
#plt.ylim(-2,2)
#plt.xlim(-.8,0)

plt.ylim(-1.2,1.2)
# plt.xlim(-.8,0)
plt.show()
#plt.savefig(FileName + "jpg", quality=100, optimize=True, dpi=600)
# plt.savefig(FileName + "jpg", quality=100, optimize=True, dpi=600)

