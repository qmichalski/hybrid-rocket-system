# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 17:36:16 2026

@author: arun2
"""
import math
import numpy as np

def property_matrix_gen(size, constant_1, constant_2, constant_3):
    base_matrix = np.zeros([size,size])
    #print(base_matrix)
    for row in range(0, size):
        if row != 0 and row != size-1:
            base_matrix[row,row - 1] = constant_3
            base_matrix[row,row]     = constant_2
            base_matrix[row,row + 1] = constant_1
        elif row == 0:
            base_matrix[row,row]     = constant_2
            base_matrix[row,row + 1] = constant_1 + constant_3
        elif row == size-1:
            base_matrix[row,row]     = constant_2
            base_matrix[row,row - 1] = constant_1 + constant_3
    #print(base_matrix)
    return base_matrix

A = np.array([[6, 1, 1],
              [4, -2, 5],
              [2, 8, 7]])
"""
A_dash = np.linalg.inv(A)
I = np.dot(A_dash,A)
print(I)


a = 1
b = 2
c = 3

N = 10 # Element number for discretization


print (property_matrix_gen(18, 9, 18, 9))
"""

T = 2000 # Exit static temperature
M = 2.56 # Exit mach number
r = 0.9 # Recovery Ratio
gamma = 1.2 # Ratio of specific heats
D = 1 # Characteristic diameter
k = 1 # Thermal conductivity
Pr = 4*gamma/(9*gamma - 5) # Prandtl Number
rho_metal = 1 # Density of metal
c_metal = 1 # Heat capacity of metal


"""
f = (0.79*(math.log(Re)) - 1.64)**(-2) # friction coefficient
Nu = (f/8)*(Re - 1000)*Pr/(1+12.7*((f/8)**(1/2))*(Pr**(2/3) - 1))
"""
Tw = 298 # Wall Temperature
T_c = 3000 # Chamber Temperature
P_c = 30*101325 # Chamber Pressure
g = 9.81 # acceleration due to gravity
D_throat = 19/1000 # Throat Diameter
visc = 4.3986/100000 # Viscosity
w = 0.6 # power for viscosity - temperature relation
C_star = 0 # Characteristic velocity
r_throat_curve = D_throat # Simplification used to make corresponding factor as unity

sigma = (((Tw*0.5/T_c)*(1 + (gamma-1)/2*M**2) + 0.5)**(0.8 - w/5))*((1 + (gamma-1)/2*M**2)**(w/5))
T_recovery = T*(1 + r*(gamma - 1)*M*M/2) # Recovery Temperature
hg = 0.026/(D_throat**0.2)*((visc**0.2)/(Pr**0.6))**0.8*(P_c*g/C_star)*((D_throat/r_throat_curve)**0.2)*sigma
alpha = k*rho_metal/c_metal




print(hg)