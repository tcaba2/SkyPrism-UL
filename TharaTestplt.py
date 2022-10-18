#!/usr/bin/env python
# coding=utf-8
import sys
import os
from importlib import *

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

src_pos_x = 2
src_pos_y = 2
map_range_x = 1
map_range_y = 1
pixels_x = 3
pixels_y = 3

ul = [[1.79028294e-10, 6.51716238e-12, 4.41879647e-12], 
[2.92180811e-12, 1.07065220e-12, 1.24647721e-12], 
[1.21292448e-12, 2.28728762e-12, 1.44649164e-12]]

dx = np.linspace(src_pos_x-map_range_x, src_pos_x+map_range_x,num=pixels_x)
dy = np.linspace(src_pos_y-map_range_y, src_pos_y+map_range_y,num=pixels_y)
dxx, dyy = np.meshgrid(dx, dy)
#dxx = dxx.reshape(pixels_x*pixels_x, 1)
#dyy = dyy.reshape(pixels_y*pixels_y, 1)

plt.pcolormesh(dxx,dyy,ul, cmap=cm.cividis)
plt.colorbar()
plt.show() 


