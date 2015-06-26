import numpy as np
from scipy.signal import fftconvolve as fftconvolve
from scipy.interpolate import interp1d
import matplotlib.pylab as plt


def gaussian(lamb,sig):
    fxn = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*((lamb-np.median(lamb))/sig)**2)
    return fxn/sum(fxn)

def sumgaus(x,lamb,nsat):
    amps = x[6:len(x)]                           #amps of satellites that can vary
    sig = x[2]                                   #sigma of central that can vary
    sigsat = sig/(nsat)                          #fixed sigma of satellites
    mult = 2.35 * sigsat/ ((lamb[-1]-lamb[0])/len(lamb))       #calculate FWHM to pixel
    shifts = ((np.linspace(0,nsat-1,nsat)) - (nsat-1)/2.)*mult #shifts of sat. in pix (one FWHM)
    fxn = gaussian(lamb,sig)                     #define center gaussian
#    plt.figure()
#    plt.plot(fxn)
    for i in range(nsat):                        #add satellites
        fxn += amps[i]*np.roll(gaussian(lamb,sigsat),int(round(shifts[i])))
#        plt.plot(amps[i]*np.roll(gaussian(lamb,sigsat),int(round(shifts[i]))))
    print x
    return fxn/sum(fxn)

def tapafunc(x,realspec0,spec,pixels,test_lam,nsat,chi_ev):
    'pull out params'
    tau = x[0]
    norm = x[1]
    sig = x[2]
    dxa = x[3]
    dxb = x[4]
    dxc = x[5]
    amps = x[6:len(x)]
    'lsf and lambda'
#    lsf = sumgaus(x,spec[:,0],nsat) 
    lsf = gaussian(spec[:,0],sig)
    lambd = dxa*pixels**2 + dxb * pixels + dxc
    'convolve tapa spec after x by tau & norm'
    model = fftconvolve(norm*(spec[:,1]**tau),lsf,'same')
    'linear interpolation'
    modelinterp = interp1d(spec[:,0],model)
    if (min(lambd) < min(spec[:,0])) or (max(lambd) > max(spec[:,0])):
        lambd = test_lam*10**3
        
    chi_ev.append(sum((realspec0 - modelinterp(lambd))**2))

#    plt.plot(lambd,realspec0,'k',lambd,modelinterp(lambd),'g')
#    print x
#    print sum((realspec0 - modelinterp(lambd))**2)

    return sum((realspec0 - modelinterp(lambd))**2)
