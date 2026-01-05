# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 10:49:22 2026

@author: arun2
"""

import CoolProp.CoolProp as CP
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')

# %% General Heat Transfer Function between any chosen fluid/solid and surface
def Q_dot_funct (rho,
                 cp,
                 c,
                 n,
                 L,
                 A,
                 g,
                 beta,
                 k,
                 T_fluid,
                 T_surface):
    
    delta_T = T_fluid - T_surface
    Ra = cp * (rho**2) * g * beta * abs(delta_T) * (L**3)
    Nu = c*Ra*n
    h = Nu * k/L
    q_dot_exchange = h * A * delta_T
    return q_dot_exchange

# %% Mass flow rate function
def mass_dot_evap(HEOS, c, a, n, L, A, g, P_tank, T_vap, T_fluid):
    
    Q_dot_liq_TO_surf = 1 
    return 0

# %% Function for rate of change of vapour volume   
def volume_dot_vap_funct(HEOS, 
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
