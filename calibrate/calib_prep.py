#!/usr/bin/python
#execfile('calib_prep.py')
import matplotlib.pylab as plt
import pyfits
from scipy.interpolate import interp1d
plt.ion()
import numpy as np


# ~ ~ ~ ~ LOAD DATA ~ ~ ~ ~ ~ #                                                           
""" spec is tapas data takes forever to load. realspec is part of the real data"""

# ~ ~ ~ ~ ~ ~ ~ ~ ~TAPAS DATA ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#Tapas Data (eventually get new spectrum with better chosen values)                      
#spec_all = np.loadtxt('/Users/ashbake/Documents/PythonProjects/GradResearch/RealData/calibrate/tapasdata.txt')                                                
istart = np.max(np.where(spec_all[:,0] > 815.))
iend = np.max(np.where(spec_all[:,0] > 835.))
spec = spec_all[iend:istart,:]
spec = spec[::-1]   #reverse array so lower wavelengths in beginning

#resize:
spec = spec[0:len(spec[:,0])/3]



# ~ ~ ~ ~ ~ ~~  REAL DATA ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
#Real Data (only one of the files)                                                      
f = pyfits.open('/Users/ashbake/Documents/PythonProjects/GradResearch/RealData/data_for_ashley2.fits')
realspecdat = f[1].data
realspec = realspecdat['FLUX'][0][ashspec]  #choose ashspec, which real spec you want

#for i in range(np.shape(realspec)[0]):
realspec = realspec/np.max(realspec)
#use dictionary eventually

errs = realspecdat['VAR'][0][ashspec]

pixels = np.arange(len(realspec))
test_lam =  (0.8171 + 0.11079*10**(-4) * (pixels) - 3.95*10**(-10)* (pixels**2))*10**3

#resize:
newsize = len(realspec)/4
realspec = realspec[0:newsize]
pixels = pixels[0:newsize]
test_lam = test_lam[0:newsize]
errs = errs[0:newsize]



# ~ ~ ~ ~ ~ ~ ~ DEFINE yy_sm ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#Define Initial Amplitude Values (to match continuum)                            
nnorm = 7
stepn = len(realspec)/nnorm
plt.figure()
xx_smalls = np.zeros(nnorm)
yy_smalls = np.zeros(nnorm)
weights = np.zeros(nnorm)
j=0
for i in range(nnorm):
    if i < (nnorm-1):
        xx = test_lam[i*stepn:(i+1)*stepn]
        yy = realspec[i*stepn:(i+1)*stepn]
    else:
        xx = test_lam[i*stepn:len(test_lam)]
        yy = realspec[i*stepn:len(test_lam)]
    xx_sm = xx[np.where((yy - np.median(yy)) > 0)] #remove large variations
    yy_sm = yy[np.where((yy - np.median(yy)) > 0)]
    if np.median(xx_sm) < 822.5 or np.median(xx_sm) > 823.5:
        xx_smalls[j] = np.median(xx_sm)
        yy_smalls[j] = np.median(yy_sm)
        weights[j] = np.std(yy_sm)
        j+=1
    plt.plot(xx_sm,yy_sm)

plt.plot(xx_smalls,yy_smalls,'go')
contin = interp1d(xx_smalls,yy_smalls,kind='cubic',bounds_error=False)
plt.plot(test_lam,contin(test_lam),'-')
