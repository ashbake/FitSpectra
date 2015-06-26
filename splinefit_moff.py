#!/usr/bin/python
#    execfile('splinefit_moff.py')
#REAL DATA VERSION
execfile('splinewrapper_moff.py')   #run code to make file
import nlopt
import time
from scipy.signal import gaussian as gaus
from scipy.signal import correlate as correlate
import functions.functions_moff as funcs #funcs.optfunc                                       


#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ NLOPT TIME ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
print funcs.bcolors.WARNING + "Starting To Do Nlopt Stuff" + funcs.bcolors.ENDC
opt = nlopt.opt(nlopt.GN_CRS2_LM,len(startvals))  #GN_CRS2_LM, MLSL
opt.set_lower_bounds(lower)
opt.set_upper_bounds(upper)
opt.set_min_objective(lambda x,grad: funcs.optfunc(x,lam_all,flux_all,N,front,end,nspec,errs))
opt.set_stopval(1e-10)
opt.set_maxeval(niter)
starttime = time.time() #time function
x = opt.optimize(startvals)
elapsedtime = time.time() - starttime
minf = opt.last_optimum_value()

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ PRINT TO TERMINAL ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

print funcs.bcolors.FAIL + '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  ' + funcs.bcolors.ENDC
print funcs.bcolors.UNDERLINE + "optimum at" + funcs.bcolors.ENDC ,x[0],x[1],x[2]
print funcs.bcolors.UNDERLINE +"minimum value = " + funcs.bcolors.ENDC, minf
print funcs.bcolors.UNDERLINE +"result code = " + funcs.bcolors.ENDC, opt.last_optimize_result()
print funcs.bcolors.UNDERLINE +"elapsed time = " + funcs.bcolors.ENDC, np.round(elapsedtime,4), ' seconds'
print funcs.bcolors.FAIL + '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ' + funcs.bcolors.ENDC

############### ~ ~ ~ ~ PLOT ~ ~ ~ ~ #######################
coeff,knots,dx,norm,tau,theta,beta,fit,small_fit = funcs.get_vals(x,lam_all,flux_all,N,front,end,nspec,numpt)
longpix, lamb, lamb2 = funcs.plt_vals(lam_all,flux_all,fit,10,knots,dx,numpt,norm,small_fit)
tck = (knots,coeff,3)


