import numpy as np
import matplotlib.pyplot as plt
import csv
FileName = "EigenValues_M800_Tres127"
Path = "D:\Documents\Codebase\Projects\TST_MKL_Eigen\TST\\"
with open(Path+FileName + '.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for index_1, row in enumerate(spamreader):
        if index_1 > 1:
            t = float(row[0])
            x, y = [], []
            for element in row[1:]:
                x.append(t)
                y.append(float(element))

            plt.plot(x, y, marker=".", linewidth=0, markersize=0.1)

plt.show()
#plt.savefig(FileName + ".jpg", quality=100, optimize=True, dpi=600)