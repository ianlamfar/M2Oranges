# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 21:55:18 2021

@author: Ian
"""

import numpy as np
import matplotlib.pyplot as plt

n = 15
z = np.pi/2
x, y = np.meshgrid(np.linspace(-2*np.pi, 2*np.pi, n), np.linspace(-2 *
                                                                  np.pi, 2*np.pi, n))

# u = 5*np.cos(x+(np.pi/2))*np.sin(y+(np.pi/2))*np.cos(z+(np.pi/2))
# v = 5*np.sin(x+(np.pi/2))*np.cos(y+(np.pi/2))*np.sin(z+(np.pi/2))
u = 5*np.cos(x+(np.pi/2))*np.sin(y+(np.pi/2))
v = 5*np.sin(x+(np.pi/2))*np.cos(y+(np.pi/2))
M = np.sqrt(u**2+v**2)

plt.figure()
a = plt.quiver(x, y, u, v, M, cmap=plt.cm.jet)
plt.colorbar(a, cmap=plt.cm.jet)
plt.show()
