#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single column radiative convective model of a dry troposphere.

The atmoshpere is well mixed and convection happens instantaneously.

Assignment 3 

1) Yes, since the variable is returned
2) The value of the albedo will be the same as the old value since the
value change inside the function has no effect on the outside value. Call by value effect?
This effect is the same if the variable is not given as an argument.
6) Changing the atmospheric temperature by the factor of 20 increases the iterations.
If the atmospheric temperature is higher the delta increases more and therefore takes longer to drop 
under the threshold. In other words, it takes longer for the system to stabilze.
7) In my opinion the maximum number of iterations which is defined by the threshold
depends on the computing capacity and the needs for precise output. With my laptop it would
be no problem to set the threshold much stricter.





"""


import numpy as np
import matplotlib.pyplot as plt


def surface_shortwave(solarIrr, albedo):
    """Calculate the sortwave radiation at the surface after albedo efects.
    S = S_0*(1-a)/4
    
    Arguments:
        solarIrr - Incoming solar radiation
        albedo - surface albedo
        
    Returns
        shortwave - absorbed shortwave radiation
    """
    shortwave = solarIrr*(1 - albedo)/4.
    return shortwave

def surface_longwave(sigma,tempAtmEm):
    longwave = sigma*tempAtmEm**4 
    return longwave

def surface_temp(incSW,incLW,sigma):
    netIncoming = incSW + incLW            # Net radiation
    tempSurf = (netIncoming/sigma)**(1./4.)# Calculate surface temperature.
    return tempSurf

def temp_due_to_convection(tempSurf,gamma,z):
    tempAtm = tempSurf - gamma*z 
    return tempAtm  

# Physical constants
solarIrr    = 1360.         # Incoming solar ratiation at the top of the atmosphere
albedo      = 0.4           # Albedo of the surface of the earth
sigma       = 5.670367E-8   # Stefan-Boltzmann constant
c2K         = 273.15        # Conversion of C to K
grav        = 9.8           # Acceleration due to gravity
cp          = 1000          # Specific heat of dry air
gamma       = grav/cp       # Dry adiabatic lapse rate

# Initial state
tempAtmEm = 50 + c2K         # Emissive temperature of the atmosphere
z = np.arange(100,10000,100) # Height levels
tempAtm = np.ones(z.shape)*tempAtmEm  # Temperature of the atmosphere
tempSurf = 50 + c2K          # Temperature at the surface

# Main program poop
delta = 1 # Change in solution. Can be anything greater than the threshold at first.
tempSurfOut = [] # variable to store the surface temperature.
count = 0 # counts number of iteration before solution is found
maxIt = 70
while abs(delta)>1e-4: # while the change in a solution
    initial = tempSurf                     # Save the temperature at the start.
    tempSurfOut.append(tempSurf - c2K)     # Add the surface temperature to the out variable.
    incSW = surface_shortwave(solarIrr, albedo) # Calculate surface shortwave.
    incLW = surface_longwave(sigma, tempAtmEm)             # Calculate surface longwave.
    tempSurf = surface_temp(incSW,incLW,sigma) # Calculate surface temperature.
    tempAtm = temp_due_to_convection(tempSurf,gamma,z)        # Calculate atmosphere temperature.
    tempAtmEm = tempAtm.mean()            # Assign the new atmosphere emmissive temperature to be the mean.
    delta = initial - tempSurf             # Find the change in the solution.
    if count > maxIt: break
    count += 1
    

# Plot the surface temperature.
plt.plot(tempSurfOut)
plt.xlabel('Iterations')
plt.ylabel('T (°C)')
plt.title('Evolution of surface temperature')
plt.show()
print("The number of iterations before a solution was found is",count)