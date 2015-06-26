#Module Defining Optimizing Function For splinefit_v#.py

import numpy as np
from scipy import interpolate
import matplotlib.pylab as plt

def optfunc(x,nu,flux,N,front,end,nspec):
    #pull out coeff                                                 
    coef = x[0:(N+4)]
    coeff = np.concatenate((coef,np.zeros(4)),axis = 1)

    #pull out knots, append                                                          
    knots_in = x[(N+4):(N+4) + N]
    knots = np.concatenate((front,knots_in,end))

    #pull out dx's
    dx = x[(N+4) + N:(N+4)+N+nspec]

    #pull out scalings
    norm = x[(N+4) + N + nspec:len(x)]
    norms = np.transpose([norm]*len(knots_in)*2)

    #defin tck, fit
    fit = np.zeros((nspec,len(nu[0])))
    for i in range(0,nspec):
        tck = (knots,coeff,3)
        fit[i,:] = interpolate.splev(nu[i] - dx[i], tck)

    return sum(sum(abs((flux/norms - fit))))



def plt_nlopt(x,nu,flux,N,front,end,nu_big,flux_big,numpt,nspec,rands):

    coef = x[0:(N+4)]
    coeff = np.concatenate((coef,np.zeros(4)),axis = 1)
    
    #define dx
    dx = x[(N+4)+N:len(x)]

    #define knots
    knots_in = x[(N+4):(N+4) + N]
    knots = np.concatenate((front,knots_in,end))

    #pull out scalings
    norm = x[(N+4) + N + nspec:len(x)]
    norms = np.transpose([norm]*len(knots_in)*2)

    #define tck & fit

    fit = np.zeros((nspec,numpt))
    for i in range(0,nspec):
        tck = (knots,coeff,3)
        fit[i,:] = interpolate.splev(nu_big[0:numpt] + rands[i] - dx[i], tck)
    for i in range(0,nspec):
        plt.figure()
        plt.plot(nu[i]-dx[i],flux[i]/norms[i][0],'go',label='data')
        plt.plot(nu_big[0:numpt] + rands[i] - dx[i],fit[i],label='NLOPT fit')
        plt.plot(nu_big[0:numpt] + rands[i] - dx[i]
                 ,flux_big[0:numpt],'r,',label='more data')
        plt.legend(loc='best')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title('Spline Fit - Nlopt')
        plt.xlim(min(knots),max(knots))
        plt.ylim(1.3,1.45)
        for j in range(0,len(knots)):
            plt.plot([knots[j],knots[j]],[1.3,1.5],'y--')
        
    return fit, dx, coeff, knots


def plt_splev(nu,flux,nu_big,flux_big,numpt,decimate):
    
    #define tck_py and fit_py
    tck_py = interpolate.splrep(nu,flux,k=3)#,t=knots_in)
    fit_py = interpolate.splev(nu_big[0:numpt],tck_py)
    
    #similar chi squared value
    print 'Python Minimum = ', sum(abs(fit_py[decimate] - flux))
    
    #plot
    
    plt.plot(nu,flux,'go',label='data')
    plt.plot(nu_big[0:numpt],fit_py,label='Python fit')
    plt.plot(nu_big[0:numpt],flux_big[0:numpt],'r,',label='more data')
    plt.legend(loc='best')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.title('Spline Fit - Python splrep')

    for i in range(0,len(tck_py[0])):
        plt.plot([tck_py[0][i],tck_py[0][i]],[1.1,1.45],'y--')


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

