# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 21:08:46 2026

@author: arun2
"""
import math
import cantera as ct
import matplotlib.pyplot as plt
gas = ct.Solution('gri30_highT.yaml')


def root_finder_funct(start, stop, tolerence, function, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat):
    if alt == False:
        initial = start
        final = stop
        mid = initial/2 + final/2
        iter_step = 0
        #print (function(initial, alt), function(final, alt))
        while abs(function(initial/2 + final/2, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat))>tolerence:
            plt.plot(iter_step, function(mid, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat),'.')
            iter_step = iter_step+1
            mid = initial/2 +final/2
            
            if function(initial, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)*function(mid, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)<0:
                final = mid
            elif function(mid, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)*function(final, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)<0:
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
        #print (function(initial, alt), function(final, alt))
        while abs(function(initial/2 + final/2, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat))>tolerence:
            plt.plot(iter_step, function(mid, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat),'.')
            iter_step = iter_step+1
            mid = initial/2 +final/2
            
            if function(initial, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)*function(mid, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)<0:
                final = mid
            elif function(mid, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)*function(final, alt, m_dot_check, P_chamber, T_chamber, propellent, R_exit, R_throat)<0:
                initial = mid
            else:
                print("Solution does not lie in the given interval!")
                break
        plt.grid()
        return mid

        

def nozzle_funct(P_estimate, alt_output = False, m_dot_check = 0, P_chamber = 30*ct.one_atm, T_chamber = 3000, propellent = 'O2:21,N2:78', R_e = 17.5/1000, R_t = 9.43/1000):
    gas.TPX = T_chamber,P_chamber,propellent
    h_chamber = gas.enthalpy_mass
    s_chamber = gas.entropy_mass
    
    gas.SP = s_chamber,P_estimate
    h_local = gas.enthalpy_mass
    a_local = gas.sound_speed
    Residual_throat = h_chamber - (h_local + (a_local**2)/2)
    v_energy = math.sqrt(2*(h_chamber-h_local))
    
    if alt_output == True:
        R_cs = R_e
    else:
        R_cs = R_t
    
    #print(R_cs)
    A_cs = math.pi*(R_cs**2)
    #m_dot_choked = gas.density_mass*A_cs*a_local
    v_continuity = m_dot_check/(gas.density_mass*A_cs)
    Residual_nozzle = v_energy - v_continuity
    
    if alt_output == False:
        return Residual_throat
    elif alt_output == True:
        return Residual_nozzle
    elif alt_output == 'All_Data':
        return gas.T, gas.sound_speed, gas.density_mass
    elif alt_output == 'Mass Flow Rate':
        return gas.density_mass*A_cs*a_local



def nozzle_solver(R_exit, R_throat, P_chamber, P_amb, T_chamber, propellent):
    P_throat = root_finder_funct(1*ct.one_atm/10, P_chamber, 0.0001, nozzle_funct, False, 0, P_chamber, T_chamber, propellent, R_exit, R_throat)
    T_throat, dummy1, dummy2 = nozzle_funct(P_throat, 'All_Data')
    m_dot_throat = nozzle_funct(P_throat, 'Mass Flow Rate')

    #plt.show()
    #print(P_throat/ct.one_atm, T_throat, m_dot_throat)

    P_exit = root_finder_funct(1*ct.one_atm/1000, P_throat*1, 0.0001, nozzle_funct, True, m_dot_throat, P_chamber, T_chamber, propellent, R_exit, R_throat)
    m_dot_exit = m_dot_throat
    T_exit, a_exit, rho_exit = nozzle_funct(P_exit, 'All_Data')
    A_exit = math.pi*(R_exit**2)
    v_exit = m_dot_exit/(rho_exit * A_exit) 
    thrust = m_dot_exit*v_exit + (P_exit - P_amb)*A_exit
    #plt.show()
    #print(P_exit/ct.one_atm, T_exit, m_dot_throat)
    #print(thrust, A_exit, rho_exit, v_exit)
    #print(m_dot_throat/(rho_exit*A_exit))
    return thrust, P_throat, T_throat, P_exit, T_exit, v_exit, a_exit, m_dot_exit

nozzle_solver(17.5/1000, 9.43/1000, 30*ct.one_atm, ct.one_atm, 3000, 'N2:0.78,O2:0.21')