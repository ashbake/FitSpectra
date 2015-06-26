# execfile('prep_data_moff.py')

"""Load and manipulate data before sending it through NLOPT"""

import matplotlib.pylab as plt
import numpy as np
import pyfits
from astropy.io import fits
import csv
from scipy import asarray as ar
from scipy.interpolate import interp1d
plt.ion()

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ LOAD DATA ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#                          

f = pyfits.open('data_for_ashley2.fits')
f2 = pyfits.open('data_for_ashley3.fits')
spec = f[1].data  # assume the first extension is a table                 
spec2 = f2[1].data
pixels = np.arange(1451)
fluxs = np.concatenate((spec['FLUX'][0],spec2['FLUX'][0]))
errs = np.concatenate((spec['VAR'][0],spec2['VAR'][0]))

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ DEFINE THINGS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#                         
numpt = 61      #ODD number of points to work with tophat/make symmetric
startpt = 510   #point in the spectrum to start out with
niter = 80000  #number of iterations for nlopt                                      
nspec = 10      #number of spectra to work with                                         
numknot = 2     #number of data points per knot point                                         

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ APPLY TO DATA  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#                        
""" make pixel array match dimensions of flux & eyeball a shift so flux"""
pix = [pixels[0:numpt]] * nspec

""" pick out numpt of flux & divide by mean """
flux = fluxs[0:0+nspec,startpt:startpt+numpt]
flux /= np.mean(flux)
refmax = sorted(flux[0])[-5]

flux = fluxs[0:0+nspec,startpt:startpt+numpt]
flux /= np.mean(flux)
refmax = sorted(flux[0])[-5]

for i in range(nspec):
    tempflux = flux[i]/np.mean(flux[i])
    flux[i] = tempflux*refmax/(sorted(tempflux)[-5])
    #plt.plot(flux[i])

plt.figure()

""" pick out numpt/ reshape errors array (exclude last 10)"""
errs = errs[0:0+nspec,startpt:startpt+numpt]/np.mean(errs)
errs[0:len(errs[0]),0:7] = errs[0:len(errs[0]),numpt-7:numpt] = np.inf


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ SHIFT SPECTRA AND FIX NORMALIZATION ~ ~ ~ ~ ~ ~ ~ ~ #                   
"""Initial Shifts"""

#output spec from wavelength calibration
tau, norm, sig, a,b,c =  3.48537391e-01,   9.88941209e-01,   1.38819216e-02, -3.95332754e-07,   1.10756971*10**-02,   8.16905624e+02
lam = a * pix[0]**2 + b * pix[0] + c #in nanometers for 0th spectrum

#shift c for all other spectra
ref_cent = lam[np.convolve(flux[0],flux[0][::-1],1).argmax()]
shifts = np.zeros(nspec)
for i in range(nspec):
    tempflux = flux[i] - np.mean(flux[i])
    tempflux /= tempflux.std()
    tempflux -= tempflux.mean()
    temp_corr = np.convolve(flux[0],tempflux[::-1],1)
    "fit peak to cubic then find max"
    maxi = temp_corr.argmax()
    pfit = np.polyfit(lam[maxi-2:maxi+2],temp_corr[maxi-2:maxi+2],2)
    testlam = np.arange(lam[maxi-2],lam[maxi+2],0.001)
    testcorr = pfit[0]*testlam**2 + pfit[1] * testlam + pfit[2]
    """find max of quadratic fit. shift in nm"""
    shifts[i] = testlam[testcorr.argmax()] - ref_cent     #pix[0][maxi] - ref_cent
    #plt.plot(lam,temp_corr,'go')
    #plt.plot(lam,temp_corr,'k-',alpha=.5)
    #plt.plot(testlam,testcorr)

###write own correlation thing?
"""Apply Shifts"""
tempshifts = np.array([shifts]*len(pix[0]))
lam_all = [lam] * nspec
lam_all = lam_all + tempshifts.T

"""Initial Taus"""
#Brute force tau finder
tau = np.zeros(nspec)
tauar = np.arange(.9,2.5,.005)
rms = np.zeros(len(tauar))
interp_flux0 = interp1d(lam_all[0],flux[0])
for i in range(nspec):
    j = 0
    for j in range(len(tauar)):
        rms[j] = np.sum(np.abs(flux[i][30:50]**tauar[j] - interp_flux0(lam_all[i][30:50])))
    tau[i] = tauar[rms.argmin()]

flux_all = np.zeros(np.shape(flux))
for i in range(nspec):
    flux_all[i,:] = flux[i]**tau[i]

"""lam_all's now aligned to within a fraction of a nanometer. Tau's roughly aligned, but need
adjustments with normalization"""

#PLOT RESULTING SPECTRA
for i in range(nspec):
    plt.plot(lam_all[i],flux_all[i])


#Coefficients that describe the first spectrum well. use to find initial taus and norms

