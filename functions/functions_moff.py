#Module Defining Optimizing Function For splinefit_v#.py

import numpy as np
from scipy import interpolate
import matplotlib.pylab as plt
from scipy.signal import fftconvolve as fftconvolve


"""Optimizing function to pass to nlopt"""
def optfunc(x,pix,flux,N,front,end,nspec,errs):
    #pull out coeff                                                                             
    coef = x[0:(N+4)]
    coeff = np.concatenate((coef,np.zeros(4)),axis = 1)

    #pull out knots, append                                                                          
    knots_in = x[(N+4):(N+4) + N]
    knots = np.concatenate((front,knots_in,end))

    #pull out dx's                                                                               
    dx = x[(N+4)+N:(N+4)+N+nspec]

    #pull out scalings                                                                                 
    norm = x[(N+4) + N + nspec:(N+4)+N + 2*nspec]

    #pull out tau's                                                                           
    tau = x[(N+4)+N+2*nspec:(N+4) + N + 3*nspec]

    #pull out thetas                                                                                    
    theta = x[(N+4)+N+3*nspec:(N+4)+N+4*nspec]

    #pull out betas                                                                                         
    beta = x[(N+4)+N+4*nspec:(N+4)+N+5*nspec]

    #define tck, fit                                                                                        
    tck = (knots,coeff,3)
    longpix = np.linspace(pix[0][0],pix[0][-1],5*(len(pix[0])))
    longfit = np.zeros((nspec,len(longpix)))
    tophat = np.zeros(len(longpix))
    tophat[len(tophat)/2 - 2:len(tophat)/2 + 3] = 1.0
    fit = np.zeros((nspec,len(pix[0])))
    index = np.arange(len(pix[0]))*5 + 2
    for bake in range(0,nspec):
        lamb = longpix + dx[bake]
#        psf = (1 + ((longpix-np.median(longpix))/theta[bake])**2)**(-1.0*beta[bake])                  
        psf = gaussian(longpix,theta[bake])
        psf = psf/np.sum(psf)
        longfit[bake,:] = norm[bake]*fftconvolve(np.abs(interpolate.splev(lamb, tck))**tau[bake],psf,'same')
        longfit[bake][np.where(np.abs(longfit[bake]) > 100)] = 100
        fit[bake,:] = np.convolve(longfit[bake]/5.0,tophat,'same')[index]

    return sum(sum(((flux - fit)**2)/errs))


"""Function to retrieve values from optimized x array"""
def get_vals(x,pix,flux,N,front,end,nspec,errs):
    #pull out coeff                                                                                     
    coef = x[0:(N+4)]
    coeff = np.concatenate((coef,np.zeros(4)),axis = 1)

    #pull out knots, append                                                                             
    knots_in = x[(N+4):(N+4) + N]
    knots = np.concatenate((front,knots_in,end))

    #pull out dx's                                                                                       
    dx = x[(N+4)+N:(N+4)+N+nspec]

    #pull out norms
    norm = x[(N+4) + N + nspec:(N+4)+N + 2*nspec]

    #pull out taus
    tau = x[(N+4)+N+2*nspec:(N+4) + N + 3*nspec]

    #pull out betas
    theta = x[(N+4)+N+3*nspec:(N+4)+N+4*nspec]

    #pull out thetas
    beta = x[(N+4)+N+4*nspec:(N+4)+N+5*nspec]

    #define tck, fit                                                                               
    tck = (knots,coeff,3)
    longpix = np.linspace(pix[0][0],pix[0][-1],5*(len(pix[0])))
    longfit = np.zeros((nspec,len(longpix)))
    tophat = np.zeros(len(longpix))
    tophat[len(tophat)/2 - 2:len(tophat)/2 + 3] = 1.0
    fit = np.zeros((nspec,len(pix[0])))
    index = np.arange(len(pix[0]))*5 + 2
    for bake in range(0,nspec):
        lamb = longpix + dx[bake]
#        psf = np.zeros(len(longpix))
#        psf[202] = 1.0
#        psf = (1 + ((longpix-np.median(longpix))/theta[bake])**2)**(-1.0*beta[bake])
#        psf = psf/np.sum(psf)
#        psf[0:120] = psf[275:404] = 0
        psf = gaussian(longpix,theta[bake])
        longfit[bake,:] = norm[bake]*fftconvolve(np.abs(interpolate.splev(lamb, tck))**tau[bake],psf,'same')
        longfit[bake][np.where(np.abs(longfit[bake]) > 100)] = 100
        fit[bake,:] = np.convolve(longfit[bake]/5.0,tophat,'same')[index]

#    truespec = interpolate.splev(lamb, tck)
    return coeff,knots,dx,norm,tau,theta,beta,longfit,fit

"""Given the values from the optimization, plot each fit for each
spectrum. Also plot the original spectrum"""
def plt_vals(lam_all,flux,fit,nspec,knots,dx,numpt,norm,small_fit):
    for i in range(0,nspec):
        plt.figure()
        longpix = np.linspace(lam_all[0][0],lam_all[0][-1],5*(len(lam_all[0])))
        lamb = longpix + dx[i]
        lamb2 = lam_all[0] + dx[i]
        plt.plot(lamb2,flux[i]/norm[i],'go',label='data')
        plt.plot(lamb2,small_fit[i]/norm[i],label='Small Fit')
        plt.plot(lamb[10:5*numpt-10],fit[i][10:5*numpt-10]/norm[i],label='NLOPT fit')
        plt.legend(loc='best')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title('Spline Fit - Nlopt Fit #' + str(i))
        plt.xlim(min(knots),max(knots))
        plt.ylim(0.0,1.1)
        for j in range(0,len(knots)):
            plt.plot([knots[j],knots[j]],[.5,1.1],'y--')        
        plt.savefig('fit_i'+str(i)+'.png')
    return longpix,lamb,lamb2

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def gaussian(lamb,sig):
    fxn = (1/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*((lamb-np.median(lamb))/sig)**2)
    return fxn/sum(fxn)
