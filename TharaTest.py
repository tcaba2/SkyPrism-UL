#!/usr/bin/env python
# coding=utf-8
import sys
import os
from importlib import *



# =================================================================================
# This script is a simplified test to calculate upper limits with systematic errors.
# By Thara Caba
# =================================================================================

from ROOT import TRolke
import math as m
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib.patches as patches

pixels_x = 3
pixels_y = 3

Non = [[7.,810.,950.],
          [475.,141.,21.],
          [15.,6.,2.]]
Noff = [[5.33333,772.333,1043.33],
          [446.333,125.333,26.667],
          [13.3333,5,1.33333]]
ex = [[6.083491122e10,2.214768808e13,2.05541865e13],
          [3.773905187e13,5.632147193e13,1.287765514e13],
          [1.181185646e13, 4.09246225e12,4.175038405e12]]

eff = 1.0
eff_err = 0.3
alpha = 0.333

tr = TRolke(0.95)

n_ulb = []
for i in range(len(Non)):
        for j in range(len(Non)):
                if Non[i][j] < Noff[i][j]:
                        Noff = np.array(Noff)
                        Noff = Noff.astype(np.int)
                        Non[i][j] = Noff[i][j]

                tr.SetGaussBkgGaussEff(Non[i][j], Noff[i][j], eff, eff_err, m.sqrt(Noff[i][j]*alpha))
                tr.SetBounding(True)
                n_ula = tr.GetUpperLimit()
                n_ulb.append(n_ula)

n_ul = []
z = int((len(n_ulb))**(0.5))
x = 0
y = z

for i in range(z):
    n = []
    for j in range (x, y):
        n.append(n_ulb[j])
    n_ul.append(n)
    x = y
    y = y + z


print(20 * "-")

ul = np.divide(n_ul, ex)

ulTrue = [[1.79037e-10, 6.51752e-12, 4.41909e-12], [2.92196e-12, 1.0701e-12, 1.24664e-12], [1.21301e-12, 2.28742e-12, 2.54521e-12]]
print(ul)
print(20 * "-")
print("This are the right results", ulTrue)
print(20 * "-")

diff = np.divide(ulTrue, ul)*100
print("This is their difference", diff)
print(20 * "-")



