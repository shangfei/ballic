#!/usr/bin/env python2
"""
This script compares the numerical solution from lane-emden.c with the analytical solution
in the case n=1.
"""
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

"""
Load the model.
"""
#model = loadtxt('lane-emden.model')
model = loadtxt('model.txt')
z   = model[:,0]
w1 = model[:,1]

# The analytic solution for n=1.
xi = arange(0.0, pi, 0.01);
theta = sin(xi)/xi

"""
Plot theta(xi)
"""
xmax = max(max(z), max(xi))*1.1
ymax = max(max(w1), max(theta))*1.1

xlim(0,xmax)
ylim(0,ymax)
xlabel(r'$\xi$')
ylabel(r'$\theta$')

plot(z, w1, '-', color='red', linewidth=1, label='Numerical')
plot(xi, theta, '--',color='blue', linewidth=1, label='Analytic')

show()

savefig('lane-emden.png')



print "Done."
