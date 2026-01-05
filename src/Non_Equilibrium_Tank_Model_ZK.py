# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 10:49:22 2026

@author: arun2
"""
# %% Imports
import CoolProp.CoolProp as CP
import math
import cantera as ct


# %% Input Data for Testing
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')

# Tank Geometry
L_tank = 1
R_tank = 1
CS_tank = math.pi*(R_tank**2)

V_tank = L_tank*CS_tank

# Initial tank data
fill_percentage = 90/100
L_vap = L_tank*(1 - fill_percentage)
L_liq = L_tank*(fill_percentage)

# Initial Tank conditions
T_vap = 298.15
T_liq = 298.15
P_tank = 56*ct.one_atm


# Heat Transfer Constants
c_dict = {
    'liquid N2O': 0.15,
    'vapour N2O': 0.15
    }
n_dict = {
    'liquid N2O': 1/3,
    'vapour N2O': 1/3
    }



# %% General Heat Transfer Function between any chosen fluid and surface
def Q_dot_funct (HEOS_from,
                 c,
                 n,
                 L,
                 A,
                 g,
                 T_from,
                 T_to):
    
    k = HEOS_from.conductivity
    beta = HEOS_from.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    delta_T = T_from - T_to
    cp = HEOS_from.cpmass
    rho = HEOS_from.rhomass
    
    Ra = cp * (rho**2) * g * beta * abs(delta_T) * (L**3)
    Nu = c*Ra*n
    h = Nu * k/L
    q_dot_exchange = h * A * delta_T
    return q_dot_exchange

# %% Mass flow rate function
def mass_dot_evap_funct (HEOS, c, a, n, L, A, g, P_tank, T_vap, T_fluid):
    
    Q_dot_liq_TO_surf = 1
    return Q_dot_liq_TO_surf

# %% Function for rate of change of vapour volume   
def volume_dot_vap_funct (HEOS, 
                         P_tank, 
                         m_vap, 
                         m_liq, 
                         m_dot_vap, 
                         m_dot_liq, 
                         V_vap, V_liq, 
                         T_vap, T_liq, 
                         T_dot_vap, 
                         T_dot_liq):
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_liq)
    partial_dP_dT_liq = HEOS.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    partial_dP_drho_liq = HEOS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
    partial_dP_dT_vap = HEOS.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    partial_dP_drho_vap = HEOS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)

    combined_term_1 = partial_dP_dT_liq * T_dot_liq + partial_dP_dT_vap * T_dot_vap
    combined_term_2 = m_dot_liq/V_liq * partial_dP_drho_liq - m_dot_vap/V_vap * partial_dP_drho_vap
    combined_term_3 = m_liq/(V_liq**2) * partial_dP_drho_liq - m_vap/(V_vap**2) * partial_dP_drho_vap
    
    dV_dt_vap = (combined_term_1 - combined_term_2)/combined_term_3
    
    return dV_dt_vap

# %% Function for rate of change of liquid volume
def volume_dot_liq_funct(V_dot_vap):
    return -V_dot_vap
