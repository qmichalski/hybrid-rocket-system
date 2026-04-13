# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 21:08:46 2026

@author: arun2
"""

import math
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
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

# nozzle_solver(17.5/1000, 9.43/1000, 30*ct.one_atm, ct.one_atm, 3000, 'N2:0.78,O2:0.21')


def root_funct(function, start, stop, tolerence):
    initial = start
    final = stop
    mid = initial/2 + final/2
    iter_step = 0
    #print (function(initial, alt), function(final, alt))
    while abs(function(initial/2 + final/2))>tolerence:
        plt.plot(iter_step, function(mid, ),'.')
        iter_step = iter_step+1
        mid = initial/2 +final/2
        
        if function(initial)*function(mid)<0:
            final = mid
        elif function(mid)*function(final)<0:
            initial = mid
        else:
            print("Solution does not lie in the given interval!")
            break
    plt.grid()
    return mid
    return None

def nozzle_sol_funct(P_ambient, P_estimate, P_chamber, T_chamber, propellent, R_throat, R_exit):
    
    A_throat = math.pi*(R_throat**2)
    A_exit = math.pi*(R_exit**2)
    
    gas.TPX = T_chamber,P_chamber,propellent
    h_chamber = gas.enthalpy_mass
    s_chamber = gas.entropy_mass
    
    #Finding Checking for throat choking
    P_choking = 0
    for P_throat_estimate in np.linspace(P_ambient/100, P_chamber, 100):
        
        gas.SP = s_chamber, P_throat_estimate - (P_ambient/100)/1000
        h_throat_estimate = gas.enthalpy_mass
        v_throat_estimate = math.sqrt(2*abs(h_chamber-h_throat_estimate))
        rho_throat_estimate = gas.density_mass
        G_throat_estimate_previous = rho_throat_estimate*v_throat_estimate
        
        gas.SP = s_chamber, P_throat_estimate
        h_throat_estimate = gas.enthalpy_mass
        v_throat_estimate = math.sqrt(2*(h_chamber-h_throat_estimate))
        rho_throat_estimate = gas.density_mass
        G_throat_estimate = rho_throat_estimate*v_throat_estimate
        T_throat = gas.T
        
        plt.plot(P_throat_estimate/ct.one_atm,G_throat_estimate,'.',color = 'blue')
        
        gas.SP = s_chamber, P_throat_estimate + (P_ambient/100)/1000
        h_throat_estimate = gas.enthalpy_mass
        v_throat_estimate = math.sqrt(2*abs(h_chamber-h_throat_estimate))
        rho_throat_estimate = gas.density_mass
        G_throat_estimate_next = rho_throat_estimate*v_throat_estimate
        check_1 = G_throat_estimate - G_throat_estimate_previous
        check_2 = G_throat_estimate_next - G_throat_estimate
        #print(check_1,check_2)
        if check_1<0 and check_2<0:
        #if G_throat_estimate >= G_throat_estimate_previous:
            #if G_throat_estimate >= G_throat_estimate_next:
            P_choking = P_throat_estimate
            break
    P_throat = P_choking
    if P_choking>P_ambient:
        #print ("Throat will choke and the pressure in the throat is ",P_choking/ct.one_atm," Atm")
        thrust_nozzle, P_throat, T_throat, P_exit, T_exit, v_exit, a_exit, m_dot_exit = nozzle_solver(R_exit, R_throat, P_chamber, ct.one_atm, T_chamber, "CO2:17.03, H2O:9.45, N2:0.5")
    else:
        #print ("Throat will NOT choke and the pressure in the throat is ",P_choking/ct.one_atm," Atm")
        # Finding flow properties when NOT choked
        gas.SP = s_chamber, P_ambient
        h_exit_NC = gas.enthalpy_mass
        v_exit_NC = math.sqrt(2*(h_chamber-h_exit_NC))
        rho_exit_NC = gas.density_mass
        G_exit_NC = rho_exit_NC*v_exit_NC
        
        P_exit = P_ambient
        
        v_exit = v_exit_NC
        
        m_dot_exit = G_exit_NC*A_exit
        T_exit = gas.T
        a_exit = gas.sound_speed
        #rho_exit = gas.density_mass
        thrust_nozzle = m_dot_exit*v_exit_NC #, T_exit 
        
        #print(m_dot_exit*v_exit_NC, T_exit, v_exit_NC)
    
    return  thrust_nozzle, P_throat, T_throat, P_exit, T_exit, v_exit, a_exit, m_dot_exit

nozzle_sol_funct(1*ct.one_atm, 100, 1.3*ct.one_atm, 298, "N2:1", 9.43/1000, 17.5/1000)