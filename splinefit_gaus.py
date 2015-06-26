#!/usr/bin/python
#    execfile('splinefit_gaus.py')
#REAL DATA VERSION
execfile('splinewrapper_gaus.py')   #run code to make file


#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ NLOPT TIME ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
print funcs.bcolors.WARNING + "Starting To Do Nlopt Stuff" + funcs.bcolors.ENDC
opt = nlopt.opt(nlopt.GN_CRS2_LM,len(startvals))  #GN_CRS2_LM
opt.set_lower_bounds(lower)
opt.set_upper_bounds(upper)
opt.set_min_objective(lambda x,grad: funcs.optfunc(x,pix,flux,N,front,end,nspec,errs))
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
coeff,knots,dxm,dxb,norm,tau,sig,fit,tck,small_fit = funcs.get_vals(x,pix,flux,N,front,end,nspec,numpt)
longpix, lamb, lamb2 = funcs.plt_vals(pix,flux,fit,nspec,knots,dxm,dxb,numpt,norm,small_fit)


xx = np.array([  1.05133385e+00,   1.07672078e+00,   2.44469715e-01,
         7.02743210e-01,   1.14840807e+00,   9.14569644e-01,
         9.58905773e-01,   1.03142130e+00,   6.27689674e-01,
         9.08271711e-01,   9.84728850e-01,   9.26257049e-01,
         9.37094108e-01,   8.01430535e-01,   6.27811801e-03,
         4.84498668e-01,   1.06516815e+00,   8.03630591e-01,
         9.91537149e-01,   2.60908367e-01,   2.26605723e-01,
         4.93101914e-03,   8.93529501e-01,   1.96272758e-01,
         9.89503240e-01,   8.09203408e-01,   1.09967851e+00,
         3.23929742e-01,   7.28539347e-01,   1.06827047e+00,
         7.73181521e-01,   8.99729825e-01,   8.47915253e-01,
         9.75769744e-01,   1.80155137e-01,   2.97711086e-01,
         5.72398618e-01,   5.42115345e-01,   1.24454007e+00,
         7.14144868e-01,   6.27780864e-01,   7.66154430e-01,
         6.32593245e-01,   7.22003408e+03,   7.22005890e+03,
         7.22008012e+03,   7.22010117e+03,   7.22012131e+03,
         7.22015327e+03,   7.22017373e+03,   7.22019753e+03,
         7.22023681e+03,   7.22024968e+03,   7.22028572e+03,
         7.22030819e+03,   7.22033226e+03,   7.22033855e+03,
         7.22038159e+03,   7.22039615e+03,   7.22041804e+03,
         7.22045497e+03,   7.22047411e+03,   7.22048790e+03,
         7.22051822e+03,   7.22054466e+03,   7.22057259e+03,
         7.22058887e+03,   7.22062493e+03,   7.22064749e+03,
         7.22068157e+03,   7.22069462e+03,   7.22072873e+03,
         7.22074927e+03,   7.22077603e+03,   7.22080719e+03,
         7.22082937e+03,   7.22083975e+03,   7.22087095e+03,
         7.22089284e+03,   7.22093331e+03,   7.22095290e+03,
         7.22097286e+03,   1.25000000e-02,   1.25000000e-02,
         1.25000000e-02,   1.25000000e-02,   1.25000000e-02,
         1.25000000e-02,   1.25000000e-02,   1.25000000e-02,
         7.21999676e+03,   7.21999793e+03,   7.22000085e+03,
         7.22000123e+03,   7.22000125e+03,   7.21999613e+03,
         7.21999686e+03,   7.21999639e+03,   1.32733984e+00,
         1.32238405e+00,   1.31003169e+00,   1.09171602e+00,
         1.28210201e+00,   1.05589641e+00,   9.96008934e-01,
         9.99782323e-01,   5.61336443e-01,   5.81172242e-01,
         5.83684054e-01,   5.72843395e-01,   5.87050955e-01,
         3.39602861e-01,   3.38497311e-01,   3.06079471e-01,
         2.27806891e-01,   5.05999890e-01,   5.27444700e-01,
         2.41990738e-01,   4.50595528e-01,   3.99065401e-01,
         3.09459978e-01,   4.21085353e-01])

