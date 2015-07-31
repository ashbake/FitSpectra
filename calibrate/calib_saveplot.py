#execfile('calib_saveplot.py')

"""Take data from Run Calib and save it and also plot things"""

plotson = 'y'       #y yes, n no


# ~ ~ ~ ~ Plot Results ~ ~ ~ ~ ~ ~ #                                              

tau = x[0]
sig = x[1]
dxa = x[2]
dxb = x[3]
dxc = x[4]
amps = x[5:5+nsat]
norms= x[5+nsat:len(x)]
    
lambd = dxa*pixels**2 + dxb * pixels + dxc
norm = interp1d(xx_smalls,norms,kind='cubic',bounds_error=False,
                fill_value=0.0)
lsf = sumgaus(x,spec[:,0],nsat,lo_inds,hi_inds)
model = fftconvolve((spec[:,1]**tau)*norm(spec[:,0]),lsf,'same')
model = fftconvolve((spec[:,1]**tau),lsf,'same')

'linear interpolation'
modelinterp = interp1d(spec[:,0],model)
finalmodel = modelinterp(lambd)
pick0s = np.where((lambd > xx_smalls[1]) & (lambd < xx_smalls[-2]))

if plotson == 'y'
    plt.figure()
    plt.plot(lambd[pick0s],realspec[pick0s],lambd[pick0s],finalmodel[pick0s])
    
    plt.figure()
    plt.plot(lambd[pick0s],(realspec[pick0s]-finalmodel[pick0s]),'g')
    plt.plot(lambd[pick0s],realspec[pick0s]/finalmodel[pick0s])
    plt.title('Residuals')
    plt.plot([lambd[0],lambd[-1]],[0,0])
    plt.figure()
    plt.plot(chi_ev,'go')
    plt.xlabel('N_iteration')
    plt.ylabel('Function Value')



# ~ Save ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
try:
    my_vars
except NameError:
    my_vars={}

my_vars['x' + str(ashspec)] = x
my_vars['res' +str(ashspec)] = realspec-finalmodel
my_vars['lambd' +str(ashspec)] = lambd
my_vars['realspec' + str(ashspec)] = realspec



#plt.figure()                                                          
plt.plot(my_vars['lambd0'],(my_vars['res0']+my_vars['res1']+my_vars['res2']+my_vars['res3']+
                            my_vars['res4']+my_vars['res5']+my_vars['res6'])/4.)
