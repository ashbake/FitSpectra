#!/usr/bin/python                                                                      
#execfile('calib_setparams.py')

# CHOOSE INITIAL PARAMETERS                                                             
nsat = 3        #number of satellite gaussians                             
nit  = 5000

# LOW START HIGH
tau = np.array([1.6, 1.82, 2.7])
sig = np.array([0.001,0.004,0.006])
dxm = np.array([-3.94e-7,-3.94e-7,-3.94e-7])
dxb = np.array([.011,.0110720461,.0111])
dxc = np.array([816.7,816.907778,817.1])
satamp1 = np.array([0.0,0.2,0.7])
satamp2 = np.array([0.0,0.2,0.7])
satamp3 = np.array([0.0,0.2,0.7])
normamps = np.array([.9*yy_smalls, 1.05*yy_smalls, 1.25*yy_smalls])


# PREPARE NLOPT INPUT ARRAYS
all_params = np.vstack((tau,sig,dxm,dxb,dxc,satamp1,satamp2,satamp3,normamps.T))
lowarr = all_params[:,0]
startarr = all_params[:,1]
higharr = all_params[:,2]

# Watch evolution of minf
chi_ev = [0]

#Set indices to make sure outside satellite gaussians are smaller
if nsat == 3:
    lo_inds = np.array([0,1,1])
    hi_inds = np.array([2,2,3])
elif nsat == 5:
    lo_inds = np.array([0,1,2,2,2])
    hi_inds = np.array([3,3,3,4,5])

