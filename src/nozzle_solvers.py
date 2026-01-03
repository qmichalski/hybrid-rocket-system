# -*- coding: utf-8 -*-
"""
Created on Sat Jan  3 18:51:25 2026

@author: Guillaume
"""
from scipy.integrate import solve_ivp
import numpy as np
import math
from scipy.optimize import root

diameter_nozzle = 19.812e-3 # m
radius_nozzle  = diameter_nozzle/2
section_nozzle = np.pi*(radius_nozzle)**2
M_design = 2.56
radius_exit = 22.4/1000
area_exit = math.pi * radius_exit ** 2
throat_area = math.pi * radius_nozzle ** 2
gamma = 1.4
Ae_at = (area_exit/throat_area)

def mach_exit_by_area(x,gamma,Ae_at):
    machRHS = (((gamma+1)/2)**(-1*(gamma+1)/(2*(gamma-1)))*((1+((gamma-1)/2)*x**2)**((gamma+1)/(2*(gamma-1))))/x) - Ae_at
    #print(machRHS)
    #print(t)
    return(machRHS)
y0 = [Ae_at]
#t_eval = np.linspace(0,tf,1000)

def solved_mach(Ae_at,machRHS,y0):
    #print(y0)
    solve = 1
    #print(machRHS)
    if (Ae_at) == machRHS:    
        solve = -1
    return (solve)
solved_mach.terminal = True
solved_mach.direction = -1



mach_exit_solver = root(mach_exit_by_area,100,args=(gamma,Ae_at))
mach_exit = mach_exit_solver['x']
print(mach_exit)
x=mach_exit
#machRHS = (((gamma+1)/2)**(-1*(gamma+1)/(2*(gamma-1)))*((1+((gamma-1)/2)*3.198**2)**((gamma+1)/(2*(gamma-1))))/3.198) - Ae_at
machRHS = (((gamma+1)/2)**(-1*(gamma+1)/(2*(gamma-1)))*((1+((gamma-1)/2)*x**2)**((gamma+1)/(2*(gamma-1))))/x) - Ae_at
print(machRHS)
"""
solve_ivp(mach_exit_by_area,[0.1,tf],y0,
                method='RK45',events=(solved_mach), 
                args=(gamma,),
                max_step=0.001)
    
mach_exit = mach_exit_solver
machRHS = ((gamma+1)/2)**(-1*(gamma+1)/(2*(gamma-1)))*((1+((gamma-1)/2)*3.198**2)**((gamma+1)/(2*(gamma-1))))/3.198
#print(machRHS)
"""