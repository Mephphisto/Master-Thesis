import numpy as np
import matplotlib.pyplot as plt
import csv
FileName = "EigenValues_M800_Tres64"
Path = "/home/jakob/Downloads/TST_MKL_Eigen/TST/cmake-build-release-intel-2019/"
with open(Path+FileName + '.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for index_1, row in enumerate(spamreader):
        if index_1 > 1:
            t = float(row[0])
            x, y = [], []
            for element in row[1:]:
                x.append(t)
                y.append(float(element))

            plt.plot(x, y, marker=".", linewidth=0, markersize=1)
        else:
            print(row)
#plt.ylim(-2,2)
plt.show()
#plt.savefig(FileName + ".jpg", quality=100, optimize=True, dpi=600)