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

# Basic Constants
g = 9.81
universal_gas_constant = 8.314

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

# Solver Inputs
time_step = 0.01

# %% General Heat Transfer Function between any chosen fluid and surface
def Heat_dot_funct (HEOS_from,
                 c,
                 n,
                 V,
                 A,
                 g,
                 T_from,
                 T_to,
                 CS = CS_tank
                 ):
    
    k = HEOS_from.conductivity
    beta = HEOS_from.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    delta_T = T_from - T_to
    cp = HEOS_from.cpmass
    rho = HEOS_from.rhomass
    L = V/CS
    
    Ra = cp * (rho**2) * g * beta * abs(delta_T) * (L**3)
    Nu = c*Ra*n
    h = Nu * k/L
    q_dot_exchange = h * A * delta_T
    return q_dot_exchange

# %% Mass flow rate due to evaporation function
def mass_dot_evap_funct (HEOS, 
                         c, 
                         n, 
                         L, 
                         A, 
                         g, 
                         P_tank, 
                         V_vap, 
                         V_liq, 
                         T_vap, 
                         T_liq, 
                         E = 693,
                         CS = CS_tank
                         ):
    
    HEOS.update(CP.PQ_INPUTS, P_tank, 1)
    h_sat_vap = HEOS.hmass
    
    HEOS.update(CP.PQ_INPUTS, P_tank, 0)
    T_surf = HEOS.T
    h_sat_liq = HEOS.hmass
    
    h_vaporization = h_sat_liq - h_sat_vap
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_liq)
    Q_dot_liq_TO_surf = Heat_dot_funct (HEOS, c["liquid N2O"], n["liquid N2O"], V_liq, A, g, T_liq, T_surf, CS)*E
    h_liq = HEOS.hmass
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
    Q_dot_surf_TO_vap = - Heat_dot_funct (HEOS, c["vapour N2O"], n["vapour N2O"], V_vap, A, g, T_vap, T_surf, CS)
    
    m_dot_evap = (Q_dot_liq_TO_surf - Q_dot_surf_TO_vap)/ (h_vaporization + h_sat_liq - h_liq)
    
    return m_dot_evap

# %% Mass flow rate due to condensation function
def mass_dot_cond_funct (HEOS,
                         P_tank, 
                         V_vap, 
                         T_vap, 
                         T_liq, 
                         E = 693,
                         CS = CS_tank,
                         R_universal = universal_gas_constant,
                         delta_t = time_step):
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
    
    Molar_mass = HEOS.molar_mass
    Z = HEOS.compressibility_factor
    
    HEOS.update(CP.QT_INPUTS, 1, T_vap)
    P_sat_vap = HEOS.P
    
    if P_tank>P_sat_vap:
        m_dot_cond = (P_tank - P_sat_vap) * V_vap * Molar_mass/ (Z * R_universal * T_vap * delta_t)
    else:
        m_dot_cond = 0
    
    return m_dot_cond

# %% Outlet Mass Flow Rate

def mass_dot_outlet_funct (HEOS, 
                           Cd, 
                           P_upstream, 
                           T_upstream, 
                           P_downstream):
    HEOS.update(CP.QT_INPUTS, 0, T_upstream)
    P_sat = HEOS.P
    
    HEOS.update(CP.PT_INPUTS, P_upstream, T_upstream)
    rho_upstream = HEOS.rhomass
    h_upstream = HEOS.hmass
    
    HEOS.update(CP.PQ_INPUTS, P_downstream, 0)
    h_downstream = HEOS.hmass
    rho_downstream = HEOS.rhomass
    
    G_SPI = Cd*math.sqrt(2*rho_upstream*(P_upstream - P_downstream))
    G_HEM = Cd*rho_downstream*math.sqrt(2*(h_upstream - h_downstream))
    K = math.sqrt((P_upstream - P_downstream)/(P_sat - P_downstream))
    
    G = (K*G_SPI + G_HEM) / (K + 1)
    
    return G
    

# %% Function for rate of change of vapour volume   
def volume_dot_vap_funct (HEOS, 
                         P_tank, 
                         m_vap, 
                         m_liq, 
                         m_dot_vap, 
                         m_dot_liq, 
                         V_vap, 
                         V_liq, 
                         T_vap, 
                         T_liq,
                         m_dot_evap, 
                         m_dot_cond,
                         m_dot_out,
                         Q_dot_in_vap,
                         Q_dot_in_liq):
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_liq)
    partial_dP_dT_liq = HEOS.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    partial_dP_drho_liq = HEOS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    cv_liq = HEOS.cvmass
    h_evap = HEOS.hmass
    h_out = HEOS.hmass
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
    partial_dP_dT_vap = HEOS.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    partial_dP_drho_vap = HEOS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    cv_vap = HEOS.cvmass
    h_cond = HEOS.hmass

    f_2 = m_dot_liq/V_liq * partial_dP_drho_liq - m_dot_vap/V_vap * partial_dP_drho_vap
    f_3 = m_liq/(V_liq**2) * partial_dP_drho_liq - m_vap/(V_vap**2) * partial_dP_drho_vap
    
    f_liq = (-m_dot_out*h_out - m_dot_evap*h_evap + m_dot_cond*h_cond + Q_dot_in_liq)/(m_liq*cv_liq)
    f_vap = (m_dot_evap*h_evap - m_dot_cond*h_cond + Q_dot_in_vap)/(m_vap*cv_vap)
    k_liq = P_tank/(m_liq*cv_liq)
    k_vap = P_tank/(m_vap*cv_vap)
    
    #f_1 = partial_dP_dT_liq * T_dot_liq + partial_dP_dT_vap * T_dot_vap #This term has been replaced with subsequent terms to make V_dot independent of T_dot
    #f_liq = (partial_dP_dT_liq/cv_liq)*(U_dot_liq/m_liq - u_liq*m_dot_liq)
    #f_vap = (partial_dP_dT_vap/cv_vap)*(U_dot_vap/m_vap - u_vap*m_dot_vap)
    #g_liq = partial_dP_drho_liq*m_dot_liq/(cv_liq*V_liq)
    #g_vap = partial_dP_drho_vap*m_dot_vap/(cv_vap*V_vap)
    #k_liq = - m_liq*partial_dP_drho_liq/(cv_liq*(V_liq**2))
    #k_vap = - m_vap*partial_dP_drho_vap/(cv_vap*(V_vap**2))
    
    dV_dt_vap = (f_liq*partial_dP_dT_liq - f_vap*partial_dP_dT_vap - f_2)/(f_3 - k_liq*partial_dP_dT_liq - k_vap*partial_dP_dT_vap)
    
    return dV_dt_vap

# %% Function for rate of change of liquid volume
def volume_dot_liq_funct(V_dot_vap):
    return -V_dot_vap
# %% Function for rate of change of fluid phase temperature
def temperature_dot_funct (HEOS, 
                           m_phase, 
                           m_dot_phase, 
                           rho_dot_phase, 
                           P_tank, 
                           T_phase, 
                           U_dot_phase
                           ):
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_phase)
    cv_phase = HEOS.cvmass
    u_phase = HEOS.umass
    partial_du_drho_phase = HEOS.first_partial_deriv(CP.iUmass, CP.iDmass, CP.iT)
    
    T_dot_phase = (1/cv_phase)*( (U_dot_phase - u_phase*m_dot_phase)/m_phase - partial_du_drho_phase*rho_dot_phase)
    
    return T_dot_phase

# %% Function for rate of change of fluid phase density
def density_dot_funct(HEOS,
                      m_phase,
                      m_dot_phase,
                      V_phase,
                      V_dot_phase):
    rho_dot_phase = m_dot_phase/V_phase - m_phase*V_dot_phase/(V_phase**2)
    return rho_dot_phase

# %% Function for Overall Liquid Mass Flow Rate
def mass_dot_LIQ_funct(m_dot_evaporation, 
                       m_dot_condensation, 
                       m_dot_exit):
    
    #m_dot_liq_overall = mass_dot_cond_funct(HEOS, P_tank, V_vap, T_vap, T_fluid, E, CS = CS_tank, R_universal = universal_gas_constant, delta_t = time_step) - mass_dot_evap_funct(HEOS, c, n, L, A, g, P_tank, V_vap, V_liq, T_vap, T_fluid, E, CS = CS_tank)
    m_dot_liq_overall = m_dot_condensation - m_dot_evaporation - m_dot_exit
    return m_dot_liq_overall

# %% Function for Overall Liquid Mass Flow Rate
def mass_dot_VAP_funct(m_dot_evaporation, 
                       m_dot_condensation):
    
    #m_dot_liq_overall = mass_dot_cond_funct(HEOS, P_tank, V_vap, T_vap, T_fluid, E, CS = CS_tank, R_universal = universal_gas_constant, delta_t = time_step) - mass_dot_evap_funct(HEOS, c, n, L, A, g, P_tank, V_vap, V_liq, T_vap, T_fluid, E, CS = CS_tank)
    m_dot_vap_overall = m_dot_evaporation - m_dot_condensation
    return m_dot_vap_overall

# %% Function for Overall Liquid ENERGY Flow Rate
def energy_dot_vap_funct(HEOS, 
                         P_tank, 
                         T_vap, 
                         T_liq,
                         m_dot_evap, 
                         m_dot_cond, 
                         Q_dot_in_vap, 
                         V_dot_vap):
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_liq)
    h_evap = HEOS.hmass
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
    h_cond = HEOS.hmass
    
    U_dot_vap = m_dot_evap*h_evap - m_dot_cond*h_cond - P_tank*V_dot_vap + Q_dot_in_vap
    
    return U_dot_vap

# %% Function for Overall Liquid ENERGY Flow Rate
def energy_dot_liq_funct(HEOS, 
                         P_tank, 
                         T_vap, 
                         T_liq, 
                         V_liq, 
                         m_dot_evap, 
                         m_dot_cond, 
                         m_dot_out, 
                         Q_dot_in_liq, 
                         V_dot_liq, 
                         threshold = 0.01):
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_liq)
    h_evap = HEOS.hmass
    
    HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
    h_cond = HEOS.hmass
    
    if V_liq < threshold*V_liq:
        HEOS.update(CP.PT_INPUTS, P_tank, T_vap)
        h_out = HEOS.hmass
    else:
        HEOS.update(CP.PT_INPUTS, P_tank, T_liq)
        h_out = HEOS.hmass
    
    U_dot_liq = - m_dot_out*h_out - m_dot_evap*h_evap + m_dot_cond*h_cond - P_tank*V_dot_liq + Q_dot_in_liq
    
    return U_dot_liq