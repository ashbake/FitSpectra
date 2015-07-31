import numpy as np
from scipy.signal import fftconvolve as fftconvolve
from scipy.interpolate import interp1d
import matplotlib.pylab as plt


def gaussian(lamb,sig):
    fxn = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*((lamb-np.median(lamb))/sig)**2)
    return fxn#/sum(fxn)

def sumgaus(x,lamb,nsat,lo_inds,hi_inds):
    amps = x[5:5+nsat]                           #amps of satellites that can vary
    sig = x[3]                                   #sigma of central that can vary
    sigsat = sig/(nsat)                          #fixed sigma of satellites
    mult = 2.35 * sigsat/ ((lamb[-1]-lamb[0])/len(lamb))    #calculate FWHM to pixel
    shifts = ((np.linspace(0,nsat-1,nsat)) - (nsat-1)/2.)*mult #shifts of sat. in pix (one FWHM)
    fxn = gaussian(lamb,sig)                     #define center gaussian
    #plt.figure()
    #plt.plot(fxn)
    for i in range(nsat):                        #add satellites
        tempamp = np.product(amps[lo_inds[i]:hi_inds[i]])
        fxn += tempamp*np.roll(gaussian(lamb,sigsat),int(round(shifts[i])))
        #plt.plot(tempamp*np.roll(gaussian(lamb,sigsat),int(round(shifts[i]))))
#        print tempamp
    return fxn/sum(fxn)



def tapafunc(x,realspec0,spec,pixels,test_lam,nsat,chi_ev,xx_smalls,lo_ind,hi_ind):
    'pull out params'
    tau = x[0]
    sig = x[1]
    dxa = x[2]
    dxb = x[3]
    dxc = x[4]
    amps = x[5:5+nsat]
    yy_smalls= x[5+nsat:len(x)]
    'lsf and lambda'
    lsf = sumgaus(x,spec[:,0],nsat,lo_ind,hi_ind) 
    lambd = dxa*pixels**2 + dxb * pixels + dxc
    'convolve tapa spec after x by tau & norm'
    norm = interp1d(xx_smalls,yy_smalls,kind='cubic',bounds_error=False,
                  fill_value=0.0)     
    model = fftconvolve((spec[:,1]**tau)*norm(spec[:,0]),lsf,'same')
    'linear interpolation'
    modelinterp = interp1d(spec[:,0],model)
    if (min(lambd) < min(spec[:,0])) or (max(lambd) > max(spec[:,0])):
        lambd = test_lam
        
    finalmodel = modelinterp(lambd)
    pick0s = np.where((lambd > xx_smalls[1]) & (lambd < xx_smalls[-2]))
    chi_ev.append(sum((realspec0[pick0s] - finalmodel[pick0s])**2))
#    plt.plot(lambd[pick0s],realspec0[pick0s],lambd[pick0s],finalmodel[pick0s])
#    print x
#    print sum((realspec0[pick0s] - finalmodel[pick0s])**2)
    return sum((realspec0[pick0s] - finalmodel[pick0s])**2)


def binres(nnorm,xin,yin):
    stepn = len(xin)/nnorm
    plt.figure()
    xx_smalls = np.zeros(nnorm-2)
    yy_smalls = np.zeros(nnorm-2)
    weights = np.zeros(nnorm-2)
    j=0
    for i in range(nnorm):
        if i < (nnorm-1):
            xx = xin[i*stepn:(i+1)*stepn]
            yy = yin[i*stepn:(i+1)*stepn]
        else:
            xx = xin[i*stepn:len(xin)]
            yy = yin[i*stepn:len(xin)]
        xx_sm = xx[np.where(np.abs(yy - np.median(yy)) < 0.005)] #remove large variations 
        yy_sm = yy[np.where(np.abs(yy - np.median(yy)) < 0.005)]
        if np.median(xx_sm) < 822.5 or np.median(xx_sm) > 823.5:
            xx_smalls[j] = np.median(xx_sm)
            yy_smalls[j] = np.median(yy_sm)
            weights[j] = np.std(yy_sm)
            j+=1
        plt.plot(xx_sm,yy_sm)
    return xx_smalls,yy_smalls


