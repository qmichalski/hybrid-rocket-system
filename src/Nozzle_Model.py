# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 21:08:46 2026

@author: arun2
"""
import math
import cantera as ct
import matplotlib.pyplot as plt
gas = ct.Solution('gri30_highT.yaml')


def root_finder_funct(start, stop, tolerence, function, alt = False, m_dot_check = 0):
    if alt == False:
        initial = start
        final = stop
        mid = initial/2 + final/2
        iter_step = 0
        #print (function(initial, alt), function(final, alt))
        while abs(function(initial/2 + final/2, alt))>tolerence:
            plt.plot(iter_step, function(mid, alt),'.')
            iter_step = iter_step+1
            mid = initial/2 +final/2
            
            if function(initial, alt)*function(mid, alt)<0:
                final = mid
            elif function(mid, alt)*function(final, alt)<0:
                initial = mid
            else:
                print("Solution does not lie in the given interval!")
                break
        plt.grid()
        return mid
    else:
        initial = start
        final = stop
        mid = initial/2 + final/2
        iter_step = 0
        #print (function(initial, alt, m_dot_check), function(final, alt, m_dot_check))
        while abs(function(initial/2 + final/2, alt, m_dot_check))>tolerence:
            plt.plot(iter_step, function(mid, alt, m_dot_check),'.')
            iter_step = iter_step+1
            mid = initial/2 +final/2
            
            if function(initial, alt, m_dot_check)*function(mid, alt, m_dot_check)<0:
                final = mid
            elif function(mid, alt, m_dot_check)*function(final, alt, m_dot_check)<0:
                initial = mid
            else:
                print("Solution does not lie in the given interval!")
                break
        plt.grid()
        return mid

        

def nozzle_funct(P_estimate, alt_output = False, m_dot_check = 0, P_chamber = 30*ct.one_atm, T_chamber = 3000, propellent = 'O2:21,N2:78'):
    gas.TPX = T_chamber,P_chamber,propellent
    h_chamber = gas.enthalpy_mass
    s_chamber = gas.entropy_mass
    
    gas.SP = s_chamber,P_estimate
    h_local = gas.enthalpy_mass
    a_local = gas.sound_speed
    Residual_throat = h_chamber - (h_local + (a_local**2)/2)
    v_energy = math.sqrt(2*(h_chamber-h_local))
    
    if alt_output == True:
        R_cs = 17.5/1000
    else:
        R_cs = 9.43/1000
    
    #print(R_cs)
    A_cs = math.pi*(R_cs**2)
    #m_dot_choked = gas.density_mass*A_cs*a_local
    v_continuity = m_dot_check/(gas.density_mass*A_cs)
    Residual_nozzle = v_energy - v_continuity
    
    if alt_output == False:
        return Residual_throat
    elif alt_output == True:
        return Residual_nozzle
    elif alt_output == 'Temperature':
        return gas.T
    elif alt_output == 'Mass Flow Rate':
        return gas.density_mass*A_cs*a_local




P_throat = root_finder_funct(1*ct.one_atm, 30*ct.one_atm, 0.0001, nozzle_funct, False)
T_throat = nozzle_funct(P_throat, 'Temperature')
m_dot_throat = nozzle_funct(P_throat, 'Mass Flow Rate')

plt.show()
print(P_throat/ct.one_atm, T_throat, m_dot_throat)

P_exit = root_finder_funct(1*ct.one_atm/1000, P_throat*1, 0.0001, nozzle_funct, True, m_dot_throat)
m_dot_exit = nozzle_funct(P_exit, 'Mass Flow Rate')
T_exit = nozzle_funct(P_exit, 'Temperature')

plt.show()
print(P_exit/ct.one_atm, T_exit, m_dot_throat)


