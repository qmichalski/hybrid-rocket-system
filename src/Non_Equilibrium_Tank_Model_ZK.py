# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 10:49:22 2026

@author: arun2
"""
# %% Imports
import CoolProp.CoolProp as CP
import math
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


# %% Input Data for Testing
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')

# Basic Constants
g = 9.81
universal_gas_constant = 8.314

# Tank Geometry
L_tank = 1
R_tank = 1
CS_tank = math.pi*(R_tank**2)
Area = CS_tank
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
                 g,
                 T_to,
                 k,
                 L,
                 A):
    
    #k = HEOS_from.conductivity()
    T_from = HEOS_from.T()
    beta = HEOS_from.isobaric_expansion_coefficient()
    delta_T = T_from - T_to
    cp = HEOS_from.cpmass()
    rho = HEOS_from.rhomass()
    
    Ra = cp * (rho**2) * g * beta * abs(delta_T) * (L**3)
    Nu = c*(Ra)**n
    h = Nu * k/L
    q_dot_from_to = h * A * delta_T
    return q_dot_from_to

# %% Mass flow rate due to evaporation function
def mass_dot_evap_funct (HEOS_sat_vap,
                         HEOS_sat_liq,
                         HEOS_liq,
                         Q_dot_liq_TO_surf,
                         Q_dot_surf_TO_vap,
                         ):
    
    h_sat_vap = HEOS_sat_vap.hmass()
    h_sat_liq = HEOS_sat_liq.hmass()
    h_vaporization = h_sat_liq - h_sat_vap
    h_liq = HEOS_liq.hmass()
    m_dot_evap = (Q_dot_liq_TO_surf - Q_dot_surf_TO_vap)/ (h_vaporization + h_sat_liq - h_liq)
    
    return m_dot_evap

# %% Mass flow rate due to condensation function
def mass_dot_cond_funct (HEOS,
                         HEOS_vap,
                         P_tank,
                         V_vap,
                         R_universal,
                         delta_t):
    
    Molar_mass = HEOS_vap.molar_mass()
    Z = HEOS_vap.compressibility_factor()
    T_V = HEOS_vap.T()
    HEOS.update(CP.QT_INPUTS, 1, T_V)
    P_sat_vap = HEOS.p()

    if P_tank>P_sat_vap:
        m_dot_cond = (P_tank - P_sat_vap) * V_vap * Molar_mass/ (Z * R_universal * T_V * delta_t)
    else:
        m_dot_cond = 0
    
    return m_dot_cond

# %% Outlet Mass Flow Rate

def mass_dot_outlet_funct (HEOS, 
                           Cd, 
                           P_upstream, 
                           T_upstream, 
                           P_downstream,
                           m_liq,
                           U_liq,
                           rho_liq
                           ):
    
    HEOS.update(CP.QT_INPUTS, 0, T_upstream)
    P_sat = HEOS.p()
    
    HEOS.update(CP.DmassUmass_INPUTS, rho_liq, U_liq/m_liq)
    rho_upstream = HEOS.rhomass()
    h_upstream = HEOS.hmass()
    
    HEOS.update(CP.PQ_INPUTS, P_downstream, 0)
    h_downstream = HEOS.hmass()
    rho_downstream = HEOS.rhomass()
    
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
                         Q_dot_in_liq,
                         U_vap,
                         U_liq):
    
    HEOS.update(CP.DmassUmass_INPUTS, m_liq/V_liq, U_liq/m_liq)
    partial_dP_dT_liq = HEOS.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    partial_dP_drho_liq = HEOS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    cv_liq = HEOS.cvmass()
    h_evap = HEOS.hmass()
    h_out = HEOS.hmass()
    
    HEOS.update(CP.DmassUmass_INPUTS, m_vap/V_vap, U_vap/m_vap)
    partial_dP_dT_vap = HEOS.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    partial_dP_drho_vap = HEOS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    cv_vap = HEOS.cvmass()
    h_cond = HEOS.hmass()

    #f_2 = m_dot_liq/V_liq * partial_dP_drho_liq - m_dot_vap/V_vap * partial_dP_drho_vap
    #f_3 = m_liq/(V_liq**2) * partial_dP_drho_liq - m_vap/(V_vap**2) * partial_dP_drho_vap
    
    #f_liq = (-m_dot_out*h_out - m_dot_evap*h_evap + m_dot_cond*h_cond + Q_dot_in_liq)/(m_liq*cv_liq)
    #f_vap = (m_dot_evap*h_evap - m_dot_cond*h_cond + Q_dot_in_vap)/(m_vap*cv_vap)
    #k_liq = P_tank/(m_liq*cv_liq)
    #k_vap = P_tank/(m_vap*cv_vap)
    
    #f_1 = partial_dP_dT_liq * T_dot_liq + partial_dP_dT_vap * T_dot_vap #This term has been replaced with subsequent terms to make V_dot independent of T_dot
    #f_liq = (partial_dP_dT_liq/cv_liq)*(U_dot_liq/m_liq - u_liq*m_dot_liq)
    #f_vap = (partial_dP_dT_vap/cv_vap)*(U_dot_vap/m_vap - u_vap*m_dot_vap)
    #g_liq = partial_dP_drho_liq*m_dot_liq/(cv_liq*V_liq)
    #g_vap = partial_dP_drho_vap*m_dot_vap/(cv_vap*V_vap)
    #k_liq = - m_liq*partial_dP_drho_liq/(cv_liq*(V_liq**2))
    #k_vap = - m_vap*partial_dP_drho_vap/(cv_vap*(V_vap**2))
    
    #dV_dt_vap = (f_liq*partial_dP_dT_liq - f_vap*partial_dP_dT_vap - f_2)/(f_3 - k_liq*partial_dP_dT_liq - k_vap*partial_dP_dT_vap)
    
    f_vap = partial_dP_dT_vap
    g_vap = partial_dP_drho_vap
    f_liq = partial_dP_dT_liq
    g_liq = partial_dP_drho_liq
    
    k_vap = (m_dot_evap*h_evap - m_dot_cond*h_cond + Q_dot_in_vap)/(m_vap*cv_vap)
    k_liq = (-m_dot_out*h_out - m_dot_evap*h_evap + m_dot_cond*h_cond + Q_dot_in_liq)/(m_liq*cv_liq)
    
    N = f_vap*k_vap - f_liq*k_liq + g_vap*m_dot_vap/V_vap + g_liq*m_dot_liq/V_liq
    D = P_tank*(f_vap/(m_vap*cv_vap) + f_liq/(m_liq*cv_liq)) + g_vap*m_vap/(V_vap**2) + g_liq*m_liq/(V_liq**2)
    
    dV_dt_vap = N/D
    
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
                           U_dot_phase,
                           U_phase,
                           rho_phase):
    #print(rho_phase,U_phase,m_phase)
    HEOS.update(CP.DmassUmass_INPUTS, rho_phase, U_phase/m_phase)
    cv_phase = HEOS.cvmass()
    u_phase = HEOS.umass()
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
                         m_dot_evap, 
                         m_dot_cond, 
                         Q_dot_in_vap, 
                         V_dot_vap
                         ):
    
    HEOS.update(CP.PQ_INPUTS, P_tank, 0)
    h_evap = HEOS.hmass()
    
    HEOS.update(CP.PQ_INPUTS, P_tank, 1)
    h_cond = HEOS.hmass()
    
    U_dot_vap = m_dot_evap*h_evap - m_dot_cond*h_cond - P_tank*V_dot_vap + Q_dot_in_vap
    
    return U_dot_vap

# %% Function for Overall Liquid ENERGY Flow Rate
def energy_dot_liq_funct(HEOS, 
                         P_tank,
                         V_liq,
                         m_dot_evap, 
                         m_dot_cond, 
                         m_dot_out, 
                         Q_dot_in_liq, 
                         V_dot_liq, 
                         threshold=0.01):

    HEOS.update(CP.PQ_INPUTS, P_tank, 0)
    h_evap = HEOS.hmass()
    HEOS.update(CP.PQ_INPUTS, P_tank, 1)
    h_cond = HEOS.hmass()

    if V_liq < threshold*V_liq:
        HEOS.update(CP.PQ_INPUTS, P_tank, 1)
        h_out = HEOS.hmass()
    else:
        HEOS.update(CP.PQ_INPUTS, P_tank, 0)
        h_out = HEOS.hmass()
    
    U_dot_liq = - m_dot_out*h_out - m_dot_evap*h_evap + m_dot_cond*h_cond - P_tank*V_dot_liq + Q_dot_in_liq

    return U_dot_liq

# %% Pressure Function

def pressure_dot_funct(HEOS,
                   m_vap,
                   m_liq,
                   m_dot_vap,
                   m_dot_liq,
                   T_vap,
                   T_liq,
                   rho_vap,
                   rho_liq,
                   U_vap,
                   U_liq,
                   T_dot_vap, 
                   T_dot_liq):

    HEOS.update(CP.DmassUmass_INPUTS, rho_liq, U_liq/m_liq)
    partial_drho_dT_liq = HEOS.first_partial_deriv(CP.iDmass, CP.iT, CP.iP)
    partial_drho_dP_liq = HEOS.first_partial_deriv(CP.iDmass, CP.iP, CP.iT)

    HEOS.update(CP.DmassUmass_INPUTS, rho_vap, U_vap/m_vap)
    partial_drho_dT_vap = HEOS.first_partial_deriv(CP.iDmass, CP.iT, CP.iP)
    partial_drho_dP_vap = HEOS.first_partial_deriv(CP.iDmass, CP.iP, CP.iT)

    N_vap = m_dot_vap/rho_vap - m_vap*partial_drho_dT_vap/(rho_vap**2)
    D_vap = m_vap*partial_drho_dP_vap/(rho_vap**2)

    N_liq = m_dot_liq/rho_liq - m_liq*partial_drho_dT_liq/(rho_liq**2)
    D_liq = m_liq*partial_drho_dP_liq/(rho_liq**2)

    P_dot = (N_vap + N_liq)/(D_vap + D_liq)

    return P_dot

# %% Modified density function

def rho_dot_funct (HEOS_vap, 
                   HEOS_liq,
                   P_tank, 
                   m_vap, 
                   m_liq, 
                   m_dot_vap, 
                   m_dot_liq, 
                   V_vap, 
                   V_liq,
                   m_dot_evap, 
                   m_dot_cond,
                   m_dot_out,
                   Q_dot_in_vap,
                   Q_dot_in_liq):

    zeta_vap = HEOS_vap.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    gamma_vap = HEOS_vap.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    CV_vap = HEOS_vap.cvmass()
    h_cond = HEOS_vap.hmass()
    alpha_vap = 1/(CV_vap*m_vap)
    beta_vap = -HEOS_vap.first_partial_deriv(CP.iUmass, CP.iDmass, CP.iT)/CV_vap
    delta_vap = -HEOS_vap.umass()*m_dot_vap/(m_vap*CV_vap)
    neta_vap = P_tank*V_vap*V_vap/m_vap

    zeta_liq = HEOS_liq.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
    gamma_liq = HEOS_liq.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
    CV_liq = HEOS_liq.cvmass()
    h_evap = HEOS_liq.hmass()
    h_out = HEOS_liq.hmass()
    alpha_liq = 1/(CV_liq*m_liq)
    beta_liq = -HEOS_liq.first_partial_deriv(CP.iUmass, CP.iDmass, CP.iT)/CV_liq
    delta_liq = -HEOS_liq.umass()*m_dot_liq/(m_liq*CV_liq)
    F_liq = - m_dot_out*h_out - m_dot_evap*h_evap + m_dot_cond*h_cond - Q_dot_in_liq
    F_vap = m_dot_evap*h_evap - m_dot_cond*h_cond + Q_dot_in_vap
    lambda_vap = F_vap + P_tank*V_vap*m_dot_vap/m_vap
    lambda_liq = F_liq + P_tank*V_liq*m_dot_liq/m_liq
    neta_liq = P_tank*V_liq*V_liq/m_liq

    A = ((V_vap*m_dot_vap/m_vap) + (V_liq*m_dot_liq/m_liq))*m_vap/(V_vap**2)
    B = ((V_liq/V_vap)**2)*(m_vap/m_liq)

    K_vap = alpha_vap*(lambda_vap + neta_vap*A) + beta_vap*A + delta_vap
    K_liq = alpha_liq*lambda_liq + delta_liq
    J_vap = -(beta_vap + alpha_vap*neta_vap)*B
    J_liq = alpha_liq*neta_liq + beta_liq

    d_rho_liq = (zeta_vap*K_vap - zeta_liq*K_liq)/(B*gamma_vap + gamma_liq - zeta_vap*J_vap + zeta_liq*J_liq)
    d_rho_vap = A - B*d_rho_liq

    d_V_liq = (V_liq/m_liq)*(m_dot_liq - V_liq*d_rho_liq)
    d_V_vap = (V_vap/m_vap)*(m_dot_vap - V_vap*d_rho_vap)

    return d_rho_vap, d_rho_liq, d_V_vap, d_V_liq



# %% Tank Model - System of ODEs
def tank_model_funct(t,z, # Main Variables
                     M_tot, time_step, HEOS, HEOS_VAP, HEOS_LIQ, HEOS_S_VAP, HEOS_S_LIQ, threshold, E, CS_tank, c_dict, n_dict, L, A, g, k_vap, k_liq, universal_gas_constant, Cd, P_exit, injector_area): # Args
    print('ITER START__________________________________________________________________________')

# Unpacking initial data____________________________________________________________________________________________________
    mass_VAP = z[0]
    mass_LIQ = z[1]
    density_VAP = z[2]
    density_LIQ = z[3]
    i_energy_VAP = z[4]
    i_energy_LIQ = z[5]
    #mass_OUT = z[6]
    print(' | m_V = ',round(mass_VAP,4),'rho_V = ',int(density_VAP),' | U_V = ',int(i_energy_VAP))
    print(' | m_L = ',round(mass_LIQ,4),'rho_L = ',int(density_LIQ),' | U_L = ',int(i_energy_LIQ))

# Intermediate parameters___________________________________________________________________________________________________
    HEOS_LIQ.update(CP.DmassUmass_INPUTS, density_LIQ, i_energy_LIQ/mass_LIQ)
    HEOS_VAP.update(CP.DmassUmass_INPUTS, density_VAP, i_energy_VAP/mass_VAP)

    pressure_tank = (HEOS_LIQ.p() + HEOS_VAP.p())/2
    print(' | Tank P = ', round(pressure_tank/ct.one_atm,2))
    
    HEOS_S_VAP.update(CP.PQ_INPUTS, pressure_tank, 0)
    HEOS_S_LIQ.update(CP.PQ_INPUTS, pressure_tank, 1)

    temperature_VAP = HEOS_VAP.T() 
    temperature_LIQ = HEOS_LIQ.T()
    temperature_SURF = HEOS_S_LIQ.T()
    print(' | T_V = ',round(temperature_VAP,2),'| T_V = ',round(temperature_LIQ,2),'| T_V = ',round(temperature_SURF,2) )
    
    volume_VAP = mass_VAP/density_VAP
    volume_LIQ = mass_LIQ/density_LIQ

    L_VAP_SURF = volume_VAP/(2*A)
    L_LIQ_SURF = L/2 - L_VAP_SURF
    #print('P_TOTAL = ',round(pressure_tank/ct.one_atm,2),'P_L = ',round(HEOS_LIQ.p()/ct.one_atm,2),'P_V',round(HEOS_VAP.p()/ct.one_atm,2),'rho_V = ',int(density_VAP),' | m_V = ',round(mass_VAP,4))

    mass_flux_out = mass_dot_outlet_funct (HEOS, 
                                              Cd, 
                                              pressure_tank, 
                                              temperature_LIQ, 
                                              P_exit,
                                              mass_LIQ,
                                              i_energy_LIQ,
                                              density_LIQ)
    dmass_OUT = injector_area * mass_flux_out
    print('Flow Rate = ', dmass_OUT*1000, 'g/s')

# Mass transfer rate for evaporation and condensation_______________________________________________________________________
    dheat_liq_to_surf = E*Heat_dot_funct (HEOS_LIQ,
                                     0.15,
                                     1/3,
                                     g,
                                     temperature_SURF,
                                     k_liq,
                                     L_LIQ_SURF,
                                     A)
    
    dheat_surf_to_vap =  Heat_dot_funct (HEOS_VAP,
                                         0.15,
                                         1/3,
                                         g,
                                         temperature_SURF,
                                         k_vap,
                                         L_VAP_SURF,
                                         A)
    print('| Q_dot L to S = ', round(dheat_liq_to_surf,2), 'J/s','| Q_dot S to V = ', round(dheat_surf_to_vap,2), 'J/s')

# Mass transfer rate for evaporation and condensation_______________________________________________________________________

    dmass_EVAP = mass_dot_evap_funct (HEOS_S_VAP,
                             HEOS_S_LIQ,
                             HEOS_LIQ,
                             dheat_liq_to_surf,
                             dheat_surf_to_vap,
                             )
    
    dmass_COND = mass_dot_cond_funct (HEOS,
                         HEOS_VAP,
                         P_tank,
                         volume_VAP,
                         universal_gas_constant,
                         time_step
                         )
    
    print('| EVAP Rate = ', round(dmass_EVAP*1000,2), 'g/s','| COND Rate = ', round(dmass_COND*1000,2), 'g/s')

# Mass transfer rate between the liquid and vapor phases____________________________________________________________________
    dmass_VAP = mass_dot_VAP_funct (dmass_EVAP, 
                                       dmass_COND
                                       )
    dmass_LIQ = mass_dot_LIQ_funct (dmass_EVAP, 
                                       dmass_COND,
                                       dmass_OUT
                                       )

# Density and Volume change rate for liquid and vapour phases__________________________________________________________________________
    drho_VAP, drho_LIQ, dvolume_VAP, dvolume_LIQ = rho_dot_funct (HEOS_VAP, 
                                                                    HEOS_LIQ,
                                                                    pressure_tank, 
                                                                    mass_VAP, 
                                                                    mass_LIQ, 
                                                                    dmass_VAP, 
                                                                    dmass_LIQ, 
                                                                    volume_VAP, 
                                                                    volume_LIQ,
                                                                    dmass_EVAP, 
                                                                    dmass_COND,
                                                                    dmass_OUT,
                                                                    dheat_surf_to_vap,
                                                                    -dheat_liq_to_surf)

# Energy transfer rate between the liquid and vapour phases_________________________________________________________________
    denergy_VAP = energy_dot_vap_funct (HEOS, 
                                           pressure_tank, 
                                           dmass_EVAP, 
                                           dmass_COND, 
                                           dheat_surf_to_vap + dheat_liq_to_surf, 
                                           dvolume_VAP
                                           )

    denergy_LIQ = energy_dot_liq_funct (HEOS, 
                                           pressure_tank,
                                           volume_LIQ,
                                           dmass_EVAP, 
                                           dmass_COND, 
                                           dmass_OUT, 
                                           -dheat_surf_to_vap-dheat_liq_to_surf,  
                                           dvolume_LIQ, 
                                           )
    #plt.plot(t,pressure_tank/ct.one_atm,'*')
    #plt.plot(t,temperature_LIQ,'.',color = 'black')
    #plt.plot(t,temperature_VAP,'.', color = 'blue')
    #plt.plot(t,dmass_EVAP,'.',color = 'black')
    #plt.plot(t,dmass_COND,'.',color = 'red')
    plt.plot(t,dvolume_VAP,'.',color = 'black')
    plt.plot(t,dvolume_LIQ,'.',color = 'red')
    #plt.plot(t,dvolume_VAP+dvolume_LIQ,'.',color = 'green')
    #plt.plot(t,dheat_liq_to_surf,'.',color = 'black')
    #plt.plot(t,dheat_surf_to_vap,'.',color = 'red')
    #plt.plot(t,mass_VAP,'.',color = 'black')
    #plt.plot(t,mass_LIQ,'.',color = 'red')
    #plt.yscale('log')
    #plt.plot(t,M_tot - (dmass_OUT*time_step + mass_VAP + mass_LIQ),'.',color = 'green')
    plt.grid()
    print(M_tot,'Kg')
    dZdt = np.zeros(len(z))
    dZdt[0] = dmass_VAP
    dZdt[1] = dmass_LIQ
    dZdt[2] = drho_VAP
    dZdt[3] = drho_LIQ
    dZdt[4] = denergy_VAP
    dZdt[5] = denergy_LIQ
    #dZdt[6] = dmass_OUT
    #print (dZdt)
    return dZdt