# execfile('splinewrapper_moff.py')
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate
execfile('prep_data_moff.py')
#lam_all and flux_all outputted

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ DEFINE INITIAL CONDITIONS, BOUNDS  ~ ~ ~ ~ ~ ~ ~ #

""" KNOT POINTS 7221.39 + 0.033756 x i 8162"""
minwav = min(lam_all[0])
maxwav = max(lam_all[0])
nknots = int(round(len(lam_all[0]))/2) #1 knot every 2 data points                                
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
hi_coef = [1.3] * (N+4)
start_coef = [1.2] * (N+4) #both should be N + 8 long at end

""" Wavelength Calibration m i + b """
mguess = (maxwav - minwav)/(numpt-1)
lo_dx = tempshifts[0] #set to match shifts found in prep_data_moff.py
hi_dx = tempshifts[0]
start_dx = tempshifts[0]

#""" Flux Calibration """
lo_scal = [.8] * nspec
hi_scal = [1.1]*  nspec
start_scal = [1.0] * nspec

#""" Tau - exp(-tau*ln(model)) """
lo_tau = [0.7] * nspec
hi_tau = [1.0] * nspec
start_tau = [1.0] * nspec

#""" Sigma of PSF """
lo_theta = [.013] * nspec  #doubles as sig when use gaussian
hi_theta = [.015] * nspec  #start at 0.014 (wavelength units nm) from wavelength calib code
start_theta = [.014] * nspec

lo_beta = [5.0] * nspec
hi_beta = [5.0] * nspec
start_beta = [5.0] * nspec

""" Concatenate """
upper = np.concatenate((hi_coef,hi_knots,hi_dx,hi_scal,hi_tau,hi_theta,hi_beta))
lower = np.concatenate((lo_coef,lo_knots,lo_dx,lo_scal,lo_tau,lo_theta,lo_beta))
#startvals = np.concatenate((start_coef,knots_in,start_dx,start_scal,start_tau,
#                            start_theta,start_beta))



