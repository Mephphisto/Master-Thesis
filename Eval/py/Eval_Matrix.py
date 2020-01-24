import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

FileName = "Matrix320_Tres64"
Path = "/home/jakob/Downloads/TST_MKL_Eigen/TST/cmake-build-release-intel-2019/"
a = []
with open(Path+FileName + '.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for index_1, row in enumerate(spamreader):
        if index_1 > 1:
            row_real = []
            for val in row[1:]:
                val = val.replace("(","")
                val = val.replace(")", "")
                row_real.append(float(val))
            a.append(row_real)
        else:
            print(row)

# 'nearest' interpolation - faithful but blocky
plt.imshow(a, interpolation='nearest')

plt.show(dpi=1200)
#plt.savefig(FileName + ".jpg", quality=100, optimize=True, dpi=600)