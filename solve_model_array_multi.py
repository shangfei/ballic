#!/usr/bin/env python2
"""
This script plots M(rho,u).
"""
from matplotlib import *
rcParams['font.family'] = 'serif'
from matplotlib.patches import *
from numpy import *
from matplotlib.pyplot import *

# This might be useful later to debug our look up table.
#isentrope = loadtxt('/home/ics/creinh/code/ic/adiabat.txt')
#data = loadtxt('modelsolve.txt',skiprows=3)
data1 = loadtxt('solve_model_array_multi.txt')
data2 = loadtxt('solve_model_array_multi.mass')

"""
Plot the models, that are within a certain range of the desired parameters.
"""
imshow(data1)
colorbar()

title('(M-Mtot)/Mtot')
xlabel('Internal energy')
ylabel('Density')

#show()
savefig('solve_model_array_multi.png')
close()

"""
Plot the mass of each model.
"""
imshow(data2)
colorbar()

title('(M-Mtot)/Mtot')
xlabel('Internal energy')
ylabel('Density')

savefig('solve_model_array_multi.mass.png')
#show()

