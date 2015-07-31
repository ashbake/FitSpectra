#!/usr/bin/python
#execfile('calib_run.py')
#Jun 22nd 2015
import os  
import nlopt
import time
from scipy.signal import gaussian as gaus
from scipy.signal import fftconvolve as fftconvolve
from  scipy.interpolate import splrep as spline
from  scipy.interpolate import splev as splev


###############################################################################
# This code optimizes wavelength calibration values (dxa,b,c), LSF, optical 
# depth (tau), continuum level (norm) to make the real spectra match the tapas
# spectrum. A multigaussian LSF is used with variable amplitudes, fixed sigma
# satellite gaussians, and central gaussian w/ variable sigma
# 
# 
#
#
##############################################################################

################################################################
# Load Data
execfile('calib_prep.py')  #only need to run first time
# Important outputs:
# spec (lamb and transmittance of tapas spec)
# testlamb,realspec,errs (real data lambda, transm., and errors
#
################################################################

#Set and Define Initial Parameters
execfile('calib_setparams.py')                                                                  
#
###############################

#Import functions
import calib_fxns
reload(calib_fxns)
from calib_fxns import gaussian,sumgaus,tapafunc,binres
###########################


############################################
#~ ~ ~ ~ ~ ~ ~ ~ Run NLOPT ~ ~ ~ ~ ~ ~ ~
opt = nlopt.opt(nlopt.GN_CRS2_LM,len(startarr)) #GN_CRS2_LM, MLSL(_LDS),GN_ESCH
opt.set_lower_bounds(lowarr)
opt.set_upper_bounds(higharr)
opt.set_min_objective(lambda x,grad: tapafunc(x,realspec,spec,pixels,test_lam,
                                              nsat,chi_ev,xx_smalls,lo_inds,hi_inds))
#opt.add_inequality_constraint(lambda x,grad: ampconstraint(x,2,1,nsat), 1e-8) #only for certain algorithms
opt.set_stopval(.0001)
opt.set_maxeval(nit)
starttime = time.time() #time function

#######################################
print 'NLOPT stuff is starting now'
x = opt.optimize(startarr)
elapsedtime = time.time() - starttime
minf = opt.last_optimum_value()

x
print minf
print elapsedtime, '  seconds'
####################################


##############################
#PLOT AND SAVE RESULTS
execfile('calib_saveplot.py')
##########################
