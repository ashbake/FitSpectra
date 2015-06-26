#!/usr/bin/python
#execfile('calib_tapas.py')
#Jun 22nd 2015
import os  
import pyfits
import nlopt
import time
import matplotlib.pylab as plt
import numpy as np
from scipy.signal import gaussian as gaus
from scipy.signal import fftconvolve as fftconvolve
from scipy.interpolate import interp1d
import calib_fxns
reload(calib_fxns)
from calib_fxns import gaussian,sumgaus,tapafunc
plt.ion()


###############################################################################
# This code optimizes wavelength calibration values (dxa,b,c), LSF, optical 
# depth (tau), continuum level (norm) to make the real spectra match the tapas
# spectrum. A multigaussian LSF is used with variable amplitudes, fixed sigma
# satellite gaussians, and central gaussian w/ variable sigma
# 
# spec:        tapas data
# realspec0:   first spectrum from the ~80 real spectra
# 
#
#
##############################################################################

nsat = 5        #number of satellite gaussians
nit  = 500

# ~ ~ ~ ~ LOAD DATA ~ ~ ~ ~ ~ #
""" spec is tapas data takes forever to load. realspec is part of the real data"""

#Tapas Data (eventually get new spectrum with better chosen values) 
#spec_all = np.loadtxt('/Users/ashbake/Documents/PythonProjects/GradResearch/RealData/calibrate/tapasdata.txt')
istart = np.max(np.where(spec_all[:,0] > 815.))
iend = np.max(np.where(spec_all[:,0] > 835.))
spec = spec_all[iend:istart,:]
spec = spec[::-1]

#Real Data (only one of the files)
f = pyfits.open('/Users/ashbake/Documents/PythonProjects/GradResearch/RealData/data_for_ashley2.fits')
realspec = f[1].data
realspec0 = realspec['FLUX'][0][0]
realspec0 = realspec0/np.max(realspec0)
errs = realspec['VAR'][0][0]

pixels = np.arange(len(realspec0))
test_lam =  0.8171 + 0.11079*10**(-4) * (pixels) - 3.95*10**(-10)* (pixels**2)

# ~ ~ ~ ~ ~ ~ ~ ~ FLATTEN realspec FLUX ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
"""split up into chunks and make median (?) of each 1"""

for i in range(21):
    if i < 20:
        xx = test_lam[i*69:(i+1)*69]
        yy = realspec0[i*69:(i+1)*69]
    else:
        xx = test_lam[i*69:len(test_lam)]
        yy = realspec0[i*69:len(test_lam)]
    xx_sm = xx[np.where((yy - np.median(yy)) > 0)] #remove large variations
    yy_sm = yy[np.where((yy - np.median(yy)) > 0)]
    a,b,c =  np.polyfit(xx_sm,yy_sm,2)
    f = a*xx**2 + b*xx + c
    if i < 20:
        realspec0[i*69:(i+1)*69] /= f
    else:
        realspec0[i*69:len(test_lam)] /=f


# ~ ~ ~ ~ ~ ~ Initial Conditions ~ ~ ~ ~ ~ ~ #
#startarr = np.array([2.5, 1.03, .01, -3.95*10**(-7),0.011079,816.9,.05,.05,.05,.05,.05]) #tau,norm,sig,dxm,dxb,dxc,amps (in fraction of 1, amp of central)
low = np.array([2.4,1.01,.008,-3.96*10**-7,.010,816.85,0,0,0,0,0])
high = np.array([2.6,1.04,.02,-3.94*10**-7,0.013,816.95,.1,.1,.1,.1,.1])
#chi_ev = [0]  #watch evolution of minf

#~ ~ ~ ~ ~ ~ ~ ~ Run NLOPT ~ ~ ~ ~ ~ ~ ~
opt = nlopt.opt(nlopt.GN_CRS2_LM,len(startarr)) #GN_CRS2_LM, MLSL(_LDS),GN_ESCH
opt.set_lower_bounds(low)
opt.set_upper_bounds(high)
opt.set_min_objective(lambda x,grad: tapafunc(x,realspec0,spec,pixels,test_lam,nsat,chi_ev))
opt.set_stopval(.0001)
opt.set_maxeval(nit)
starttime = time.time() #time function                                                        
print 'NLOPT stuff is starting now'
x = opt.optimize(startarr)
elapsedtime = time.time() - starttime
minf = opt.last_optimum_value()

print x
print minf
print elapsedtime, '  seconds'

# ~ ~ ~ ~ Plot Results ~ ~ ~ ~ ~ ~ #
lambd = x[3] * pixels**2 + x[4] * pixels + x[5]
plt.figure()
lsf = sumgaus(x,spec[:,0],nsat)

#lsf = gaussian(spec[:,0],sig)
tau = x[0]
norm = x[1]
sig = x[2]
dxa = x[3]
dxb = x[4]
dxc = x[5]
amps = x[6:len(x)]
model = fftconvolve(norm*(spec[:,1]**tau),lsf,'same')
modelinterp = interp1d(spec[:,0],model)
plt.plot(lambd,realspec0,'k',lambd,modelinterp(lambd),'g')

plt.figure()
plt.plot(lambd,realspec0-modelinterp(lambd),'g')

plt.figure()
plt.plot(chi_ev)
plt.xlabel('N_iteration')
plt.ylabel('Function Value')


#best one
#array([  3.48537391e-01,   9.88941209e-01,   1.38819216e-02,
#        -3.95332754e-07,   1.10756971e-02,   8.16905624e+02])


