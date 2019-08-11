import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    dimensiony = int(sys.argv[1])
    name = sys.argv[2]
    dimensionx = (len(sys.argv) - 3) / dimensiony
    a = np.zeros((dimensionx, dimensiony))
    i = 3;
    while (i < len(sys.argv)):
        for j in range(dimensiony):
            a[(i - 3) / dimensiony][j] = float(sys.argv[i])
            i += 1
    plt.style.use('ggplot')
    if (dimensiony == 3):
        # plt.figure(figsize=(12,8),dpi=200)
        fig = plt.figure(figsize=(8, 6), dpi=100)
        ax = Axes3D(fig)
        ax.scatter(a[:, 0], a[:, 1], a[:, 2], c="b")
        plt.xlabel('x')
        plt.ylabel("y")
        plt.title(name)
        ax.view_init(32, 45)
        plt.show()
        plt.close()
    elif (dimensiony == 2):
        plt.figure(figsize=(12, 8), dpi=80)
        plt.plot(a[:, 0], a[:, 1], "*b")
        plt.xlabel('x')
        plt.ylabel("y")
        plt.title(name)
        plt.show()
    else:
        d = np.arange(1, dimensiony + 1)
        plt.figure(figsize=(12, 8), dpi=80)
        for i in range(dimensionx):
            plt.plot(d, a[i, :], "b")
        plt.xlabel('Values')
        plt.ylabel("Objctive")
        plt.xticks(range(1, dimensiony + 1))
        plt.title(name)
        plt.show()
