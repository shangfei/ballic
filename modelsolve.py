#!/usr/bin/env python2
"""
This script produces a u(rho) plot for a planetary model.
"""
#import matplotlib as mpl; mpl.rcParams['font.family'] = 'serif'
from matplotlib import *
rcParams['font.family'] = 'serif'
from numpy import *
from matplotlib.pyplot import *
from sys import *

"""
if len(argv) != 2:
		print "Usage: modelic.py model"
		exit(1)

modelfile = argv[1]
model = loadtxt(modelfile)
"""

model = loadtxt('ballic.model')
output = loadtxt('modelsolve.txt')

"""
Load the ballic model.
"""
R1   = model[:,0]
rho1 = model[:,1]
M1   = model[:,2]
u1   = model[:,3]

"""
Load the output from modelsolve.c.
"""
R2   = output[:,0]
rho2 = output[:,1]
M2   = output[:,2]
u2   = output[:,3]

"""
Plot M(R)
"""
xmax = max(max(R1), max(R2))*1.1
ymax = max(max(M1), max(M2))*1.1

xlim(0,xmax)
ylim(0,ymax)
xlabel('Radius')
ylabel('Mass')

plot(R1,M1,'-',color='red',linewidth=1,label='Model 1')
plot(R2,M2,'--',color='blue',linewidth=1,label='Model 2')
scatter(R2,M2,s=25,color='green')
savefig('modelsolve.M.png')
close()

"""
Plot rho(R)
"""
xmax = max(max(R1), max(R2))*1.1
ymax = max(max(rho1), max(rho2))*1.1

xlim(0,xmax)
ylim(0,ymax)
xlabel('Radius')
ylabel('Density')

plot(R1,rho1,'-',color='red',linewidth=1,label='ballic.model')
plot(R2,rho2,'--',color='blue',linewidth=1,label='Int. from surf.')
scatter(R2,rho2,s=10,color='green')
legend(loc='best')
show()
savefig('modelsolve.rho.png')
close()

"""
Plot u(R)
"""
xmax = max(max(R1), max(R2))*1.1
ymax = max(max(u1), max(u2))*1.1

xlim(0,xmax)
ylim(0,ymax)
xlabel('Radius')
ylabel('Internal energy')

plot(R1,u1,'-',color='red',linewidth=1,label='Model 1')
plot(R2,u2,'--',color='blue',linewidth=1,label='Model 2')

savefig('modelsolve.u.png')
close()




print "Done."
