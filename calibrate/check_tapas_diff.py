#execfile('check_tapas_diff.py')
import numpy as np
import matplotlib.pylab as plt
plt.ion()

spec1 = np.loadtxt('/Users/ashbake/Documents/PythonProjects/GradResearch/RealData/calibrate/tapas_dataDec20.txt')
spec2 = np.loadtxt('/Users/ashbake/Documents/PythonProjects/GradResearch/RealData/calibrate/tapas_dataDec20.txt')

lam = spec1[:,0]
y1 = spec1[:,1]
y2 = spec2[:,1]

plt.plot(lam,y2-y1)
