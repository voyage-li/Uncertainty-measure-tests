# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os

fig = plt.figure()
ax = Axes3D(fig)
X = np.arange(-4, 4, 0.25)
Y = np.arange(-4, 4, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = X*Y

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')

plt.draw()
try:
    plt.savefig('./temp/3D.jpg')
except:
    os.system('mkdir temp')
    plt.savefig('./temp/3D.jpg')

plt.close()
