#!/usr/bin/env python
# coding=utf-8
import sys
import os
from importlib import *

# ===========================================================================
# This script is mainly to calculate upper limits with systematic errors.
# By Thara Caba
# ===========================================================================

from ROOT import TRolke
import math as m
import numpy as np
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib.patches as patches

################## Loading SkyPrism #############

sys.path.append(os.path.expandvars("/"))

import Configuration
reload(Configuration)
from Configuration import *

import General_functions
reload(General_functions)
from General_functions import *

import SkyMap
reload(SkyMap)
from SkyMap import *

try:
  Conf = Configuration(sys.argv[1])
except:
  print "You did not specify a config file: Trying default filename 'skyprism.rc'."
  Conf = Configuration()

Projectname = Conf.Get("Projectname")
ELow  = map(float, Conf.Get('Energies').split(','))[:-1]
pixels_x = int(Conf.Get("ULMapBinsRa"))
pixels_y = int(Conf.Get("ULMapBinsDec"))

################## Loading sky maps #############

for Bin, E in enumerate(ELow):
        OnMap = SkyMap()
        OnMap.readFITS(Conf.Get("Outfilepath") + Projectname + "_OnMap" + ".fits", Bin+1)
        Non = OnMap.Image
        HadMap = SkyMap()
        HadMap.readFITS(Conf.Get("Outfilepath") + Projectname + "_BkgMap" + ".fits", Bin+1)
        Noff = HadMap.Image

        ExpMap = SkyMap()
        ExpMap.readFITS(Conf.Get("Outfilepath")+ Projectname + "_Exposure_Eest" + ".fits", Bin+1)
        ex = ExpMap.Image

################## Applying the Rolke Method #############

tr = TRolke(0.95)

eff = 1.0
eff_err = 0.3
alpha = 0.333

Non = np.array(Non)

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

ul = np.divide(n_ul, ex)

ul = ul.reshape(pixels_x*pixels_x, 1)

results_list = "./upperlimits_CrabNebula_Ebin0.txt"
results = np.loadtxt(results_list)
ts = results[:,6]

for i in range (len(ul)):
        if ts[i] > 5:
                ul[i] = 0

np.savetxt('ulTest.txt', ul, fmt='%.16e')

print("Well done! Time to plot")

################ plotting vlues ##################

# manual recalculating of bin middle to show fluxpoint as middle of a bin
ra = results[:,0]*0.066667
dec = results[:,1]

dx = np.unique(ra) #ra
dy = np.unique(dec) #dec

#  calculating middle of each bin for the map
diff_x = np.diff(dx)[0]
range_x = np.ones(len(dx) + 1)
range_x[:-1] = dx
range_x[-1] = dx[-1] + diff_x
range_x -= diff_x / 2.
range_strange_x = range_x[::-1]

diff_y = np.diff(dy)[0]
range_y = np.ones(len(dy) + 1)
range_y[:-1] = dy
range_y[-1] = dy[-1] + diff_y
range_y -= diff_y / 2.

# plot vlues
plt.hist2d(ra, dec, weights=ul[:,0], cmap=cm.YlGnBu_r, bins=(range_x, range_y), label="Flux Upper Limit", rasterized=True)
plt.xlim(max(ra), min(ra))
plt.title("Flux Upper Limits Skymap")
plt.xlabel('RA [h]')#, fontsize=labelsize)
plt.ylabel('Dec [deg]')#, fontsize=labelsize)
red_triangle = mlines.Line2D([], [], color='red', marker='^', linestyle='None', markersize=10, label='red triangle')
plt.legend(handles=[red_triangle])
plt.legend([red_triangle], ['Crab Nebula'], loc='upper left', framealpha=0.5)
red_triangle = plt.plot(5.5755, 22.0145, '^', c='red')
cbar = plt.colorbar()
cbar.set_label("Flux Upper Limits / cm$^{-2}$s$^{-1}$s", rotation=270, labelpad=28)
plt.savefig("Flux2Test")
plt.close()

plt.hist(ul, histtype='stepfilled', color='mistyrose', ec='red', bins=30)
plt.title("UL value distribution")
plt.xlabel('Flux Upper Limit')#, fontsize=labelsize)
plt.ylabel('N')#, fontsize=labelsize)
plt.savefig("fluxULDistributionTest")
plt.close()

print("Well done! Check your plots")

################ Verifying results ##################

diff = np.divide(results[:,2], ul[:,0])*100
print("This is their difference", diff)

plt.hist(results[:,2], histtype='stepfilled', color='lightblue', ec='red', bins=30, alpha=0.5)
plt.hist(ul, histtype='stepfilled', color='mistyrose', ec='blue', bins=30)
plt.title("UL value distribution")
plt.xlabel('Flux Upper Limit')#, fontsize=labelsize)
plt.ylabel('N')#, fontsize=labelsize)
plt.savefig("Comparison")
plt.close()

