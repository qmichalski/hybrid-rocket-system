# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 10:49:22 2026

@author: arun2
"""

import CoolProp.CoolProp as CP


def Q_dot_funct(rho,cp,c,n,L,A,g,beta,k,T_fluid, T_surface):
    delta_T = T_fluid - T_surface
    Ra = cp * (rho**2) * g * beta * abs(delta_T) * (L**3)
    Nu = c*Ra*n
    h = Nu * k/L
    q_dot_exchange = h * A * delta_T
    return q_dot_exchange

def mass_dot_funct():
    print ('Sample')