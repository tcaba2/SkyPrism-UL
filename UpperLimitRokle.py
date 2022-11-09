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

OKGREEN   = '\033[92m'
ENDC	  = '\033[0m'

try:
  Conf = Configuration(sys.argv[1])
except:
  print OKGREEN + "You did not specify a config file: Trying default filename 'skyprism.rc'." + ENDC
  Conf = Configuration()

Projectname = Conf.Get("Projectname")
ELow  = map(float, Conf.Get('Energies').split(','))[:-1]
pixels_x = int(Conf.Get("ULMapBinsRa"))
pixels_y = int(Conf.Get("ULMapBinsDec"))
SrcRa  = float(Conf.Get('SrcRa'))
SrcDec = float(Conf.Get('SrcDec'))
BgMethod = Conf.Get("BgMethod")

try:
  CircRegionsRa   = map(float, Conf.Get('CircRegionsRa',  optional=True).split(','))
  CircRegionsDec  = map(float, Conf.Get('CircRegionsDec', optional=True).split(','))
  CircRegionsRad  = map(float, Conf.Get('CircRegionsRad', optional=True).split(','))
except:
  CircRegionsRa   = []
  CircRegionsDec  = []
  CircRegionsRad  = []
try:
  LineRa1   = map(float, Conf.Get('LineRa1',  optional=True).split(','))
  LineDec1  = map(float, Conf.Get('LineDec1', optional=True).split(','))
  LineRa2   = map(float, Conf.Get('LineRa2',  optional=True).split(','))
  LineDec2  = map(float, Conf.Get('LineDec2', optional=True).split(','))
  LineDist  = map(float, Conf.Get('LineDist', optional=True).split(','))
except:
  LineRa1   = []
  LineDec1  = []
  LineRa2   = []
  LineDec2  = []
  LineDist  = []


################## Loading sky maps #############
print OKGREEN + "Loading sky maps" + ENDC

for Bin, E in enumerate(ELow):
        OnMap = SkyMap()
        OnMap.readFITS(Conf.Get("Outfilepath") + Projectname + "_OnMap" + ".fits", Bin+1)
        Non = OnMap.Image.T

        HadMap = SkyMap()
        HadMap.readFITS(Conf.Get("Outfilepath") + Projectname + "_BkgMap" + ".fits", Bin+1)

        ExpMap = SkyMap()
        ExpMap.readFITS(Conf.Get("Outfilepath")+ Projectname + "_Exposure_Eest" + ".fits", Bin+1)
        ex = ExpMap.Image.T

        RaDec = np.meshgrid(OnMap.GetRA_axis().GetCenter(),OnMap.GetDec_axis().GetCenter())
        Ra = RaDec[0]
        Dec = RaDec[1]
        X = (SrcRa-Ra)*m.cos(m.radians(SrcDec))
        Y = Dec-SrcDec

        #Normalizing the Off map
        ExMap = 0.0*OnMap.Image
        for i in range(len(CircRegionsRad)):
                dRA   = (Ra - CircRegionsRa[i]) * m.cos(m.radians(Dec))
                dDec  = Dec - CircRegionsDec[i]
                dist2 = (dRA**2 + dDec**2)
                ExMap += 1.0 * (dist2 < CircRegionsRad[i]**2)
        for i in range(len(LineDist)):
                x0 = (Ra           - SrcRa) * m.cos(m.radians(Dec))
                y0 =  Dec          - SrcDec
                x1 = (LineRa1 [i] - SrcRa) * m.cos(m.radians(Dec))
                y1 =  LineDec1[i] - SrcDec
                x2 = (LineRa2 [i] - SrcRa) * m.cos(m.radians(Dec))
                y2 =  LineDec2[i] - SrcDec
                dist2 = pow( (y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1, 2) / ( pow(y2-y1, 2) + pow(x2-x1, 2) )
                ExMap += 1.0 * ( dist2 < LineDist[i]**2 )
        
        if BgMethod == "WobbleMap":
                radius = 0.4
                NormMask = np.where( X**2 + Y**2 > radius**2 )

        elif BgMethod == "ExclMap" or BgMethod == "OffMap" or BgMethod == "ExclMapFit" :
                NormMask = np.where( 1.0 - ExMap.transpose() > 0 )
        else:
             	radius = 0.0
                NormMask = np.where( X**2 + Y**2 > radius**2 )

        ExpMask = ExpMap.Image > ExpMap.Image.max()/20.

        MedianRatio    = np.mean(OnMap.Image[NormMask])/np.mean(HadMap.Image[NormMask])
        N_HadMap = HadMap.Scale(MedianRatio)

        Noff = N_HadMap.Image.T

################## Applying the Rolke Method #############
print OKGREEN + "Calculating upper limits" + ENDC

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

#Making n_ul a matrix
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

results_list = "./upperlimits_%s_*.txt" % (Projectname)
results = np.loadtxt(results_list)
ts = results[:,6]

#Upper limits only below 5 sigma

for i in range (len(ul)):
        if ts[i] > 5:
                ul[i] = 0

np.savetxt('ulTest.txt', ul, fmt='%.16e')

print OKGREEN + "Well done time to plot!" + ENDC

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

#This part should be changed for each map
red_triangle = mlines.Line2D([], [], color='red', marker='^', linestyle='None', markersize=10, label='red triangle')
plt.legend(handles=[red_triangle])
plt.legend([red_triangle], ['Crab Nebula'], loc='upper left', framealpha=0.5)
red_triangle = plt.plot(5.5755, 22.0145, '^', c='red')

#This part doesn't need change for each map
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

print OKGREEN + "Well done! Check your plots" + ENDC

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

plt.hist(diff, histtype='stepfilled', color='lightblue', ec='red', bins=30, range=(0, 2000))
plt.savefig("diff")
plt.close()






