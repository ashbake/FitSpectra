#execfile('splinewrapper_gaus.py')
import nlopt
import matplotlib.pylab as plt
import numpy as np
import pyfits
from scipy.optimize import curve_fit
from scipy import asarray as ar
from astropy.io import fits
import csv
from scipy.interpolate import interp1d
from scipy import interpolate
import time
import functions.functions_gaus as funcs #funcs.optfunc
from scipy.signal import gaussian as gaus
plt.ion()

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ LOAD DATA ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#

f = pyfits.open('data_for_ashley2.fits')
spec = f[1].data  # assume the first extension is a table
pixels = np.arange(1451)
fluxs = spec['FLUX'][0]
errs = spec['VAR'][0]

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ DEFINE THINGS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
#for i in range(8):
#    plt.plot(pix[i],flux[i])

numpt = 81      #ODD number of points to work with from the spectrum
niter = 50000  #number of iterations for nlopt
nspec = 8      #number of spectra to work with
numknot = 2     #number of data points per knot point

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ APPLY TO DATA  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
""" make pixel array match dimensions of flux & eyeball a shift so flux"""
pix = [pixels[0:numpt]]*nspec
""" pick out numpt of flux & divide by mean """
flux = fluxs[2:nspec+2,510:510+numpt]
meanflux = np.mean(flux)
flux = flux/meanflux
""" pick out numpt/ reshape errors array (exclude last 10)"""
errs = errs[2:2+nspec,0:numpt]/np.mean(errs)
errs[0:len(errs[0]),0:7] = errs[0:len(errs[0]),numpt-7:numpt] = np.inf

for i in range(nspec):
    plt.plot(pix[i],flux[i])
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ DEFINE INITIAL CONDITIONS, BOUNDS  ~ ~ ~ ~ ~ ~ ~ #

""" KNOT POINTS 7221.39 + 0.033756 x i 8162"""
minwav = 7220.0
maxwav = 7221.0
nknots = int(round(len(pix[0]))/2) #1 knot every 2 data points                                
gridsize = (maxwav - minwav)/np.double(nknots)
knots_in = np.arange(minwav+gridsize,maxwav-gridsize,gridsize)
N = len(knots_in) #doesnt include 4 repeats at endpoint and 4 at front                           
front = np.array([0.0]*4) + minwav
end = np.array([0.0]*4) + maxwav
knots = np.concatenate((front,knots_in,end))
lo_knots = knots_in - gridsize/2.
hi_knots = knots_in + gridsize/2.

""" Coefficients """
lo_coef = [0.0] * (N+4) # each segment has N+4 (last 4 zero)
hi_coef = [1.5] * (N+4)
start_coef = [1.1] * (N+4) #both should be N + 8 long at end

""" Wavelength Calibration m i + b """
mguess = (maxwav - minwav)/(numpt-1)
lo_m = [mguess]*nspec #set the first spectrums m1. others should be near
hi_m = [mguess]*nspec
start_m = [mguess] * nspec

#if nspec > 1:
#    hi_m = np.concatenate((lo_m,[lo_m[0] + 0.005] * (nspec - 1)))
#    lo_m = np.concatenate((lo_m,[lo_m[0] - 0.005] * (nspec - 1)))

lo_b = [minwav - gridsize]*nspec #[minwav]
hi_b = [minwav + gridsize]*nspec  #[minwav]
start_b = [minwav] * nspec

#if nspec > 1:
#    hi_b = np.concatenate((lo_b,[lo_b[0] + 2*gridsize] * (nspec - 1)))
#    lo_b = np.concatenate((lo_b,[lo_b[0] - 2*gridsize] * (nspec - 1)))

""" Flux Calibration """
lo_scal = [0.5] * nspec
hi_scal = [1.7]*  nspec
start_scal = [1.0] * nspec

""" Tau - exp(-tau*ln(model)) """
lo_tau = [.1] * nspec
hi_tau = [1.0] * nspec
start_tau = [1.0] * nspec

#""" Sigma of PSF """
lo_sig = [.2] * nspec
hi_sig = [.4] * nspec
start_sig = [.3] * nspec

""" Concatenate """
upper = np.concatenate((hi_coef,hi_knots,hi_m,hi_b,hi_scal,hi_tau,hi_sig))
lower = np.concatenate((lo_coef,lo_knots,lo_m,lo_b,lo_scal,lo_tau,lo_sig))
startvals = np.concatenate((start_coef,knots_in,start_m,start_b,start_scal,start_tau,
                            start_sig))


#numpt 81, nspec 8, startpix 510, min val = 6.5
#edges seem shifted maybe b/c of edge effects. tau's good but sigmas don't make sense
xnice = np.array([  8.92962771e-01,   1.12527557e+00,   2.81136691e-01,
         7.41063190e-01,   9.93338134e-01,   9.50031629e-01,
         9.30675538e-01,   1.01342359e+00,   5.87188304e-01,
         9.81999420e-01,   9.25523149e-01,   9.39099949e-01,
         9.13114057e-01,   7.68279241e-01,   6.56818764e-04,
         6.24413060e-01,   9.83678848e-01,   8.30313603e-01,
         9.93007273e-01,   3.16639571e-01,   2.66804762e-01,
         5.29900161e-04,   9.27656435e-01,   1.36435158e-01,
         9.76499056e-01,   8.02936226e-01,   1.08638909e+00,
         2.74286799e-01,   8.59249751e-01,   9.42029008e-01,
         8.64974045e-01,   8.32146357e-01,   9.02545326e-01,
         8.69947266e-01,   1.23550883e-02,   4.96656354e-01,
         6.31018731e-01,   3.40353856e-01,   1.49955436e+00,
         2.09645897e-01,   5.11221825e-01,   6.66060193e-01,
         4.33291961e-01,   7.22003271e+03,   7.22005977e+03,
         7.22007967e+03,   7.22010468e+03,   7.22012280e+03,
         7.22016109e+03,   7.22018608e+03,   7.22019919e+03,
         7.22023555e+03,   7.22025403e+03,   7.22028369e+03,
         7.22031242e+03,   7.22033296e+03,   7.22033844e+03,
         7.22038187e+03,   7.22039572e+03,   7.22041271e+03,
         7.22044617e+03,   7.22047995e+03,   7.22048784e+03,
         7.22051251e+03,   7.22053751e+03,   7.22057069e+03,
         7.22058816e+03,   7.22063740e+03,   7.22065511e+03,
         7.22068086e+03,   7.22069311e+03,   7.22072719e+03,
         7.22074917e+03,   7.22077935e+03,   7.22081249e+03,
         7.22083749e+03,   7.22084620e+03,   7.22086665e+03,
         7.22088760e+03,   7.22092922e+03,   7.22095160e+03,
         7.22097063e+03,   1.25000000e-02,   1.25000000e-02,
         1.25000000e-02,   1.25000000e-02,   1.25000000e-02,
         1.25000000e-02,   1.25000000e-02,   1.25000000e-02,
         7.21999636e+03,   7.21999828e+03,   7.22000023e+03,
         7.22000010e+03,   7.22000035e+03,   7.21999621e+03,
         7.21999623e+03,   7.21999593e+03,   1.34355379e+00,
         1.32792818e+00,   1.32071229e+00,   1.10109534e+00,
         1.29356696e+00,   1.06471891e+00,   9.98873566e-01,
         1.00278061e+00,   6.09678842e-01,   6.11616592e-01,
         6.18598176e-01,   6.21963570e-01,   6.18284980e-01,
         3.68061227e-01,   3.49173723e-01,   3.18653728e-01,
         1.19795633e-01,   1.18705705e-01,   2.42183117e-01,
         2.28621073e-01,   2.27706935e-01,   1.36789897e-01,
         1.58700568e-01,   3.41365051e-01])
