# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 22:08:42 2026

@author: arun2
"""
import CoolProp.CoolProp as CP
import math
import cantera as ct
import numpy as np
import Non_Equilibrium_Tank_Model_ZK as ZK
from scipy.integrate import solve_ivp

# %% Tank Model - System of ODEs
def tank_model_funct(t,z, # Main Variables
                     HEOS, threshold, E, CS_tank, c_dict, n_dict, L, A, g, universal_gas_constant): # Args

# Unpacking initial data____________________________________________________________________________________________________
    pressure_tank = z[0]
    temperature_VAP = z[1]
    temperature_LIQ = z[2]
    mass_VAP = z[3]
    mass_LIQ = z[4]
    #mass_EVAP = z[5]
    #mass_COND = z[6]
    #mass_OUT = z[7]
    #rho_VAP = z[8]
    #rho_LIQ = z[9]
    volume_VAP = z[5]
    volume_LIQ = z[6]
    #i_energy_VAP = z[12]
    #i_energy_LIQ = z[13]
    
# Intermediate parameters___________________________________________________________________________________________________
    dmass_OUT = 0
    HEOS.update(CP.PQ_INPUTS, pressure_tank, 0)
    temperature_surf = HEOS.T()
    print(temperature_surf,temperature_LIQ)
    
# Mass transfer rate for evaporation and condensation_______________________________________________________________________
    dheat_in_liq = ZK.Heat_dot_funct (HEOS, 
                                      c_dict["liquid N2O"], 
                                      n_dict["liquid N2O"], 
                                      volume_LIQ, 
                                      A, 
                                      g, 
                                      temperature_LIQ, 
                                      temperature_surf, 
                                      CS_tank)*E
    
    dheat_in_vap = ZK.Heat_dot_funct (HEOS, 
                                      c_dict["vapour N2O"], 
                                      n_dict["vapour N2O"], 
                                      volume_VAP, 
                                      A, 
                                      g, 
                                      temperature_VAP, 
                                      temperature_surf, 
                                      CS_tank)

# Mass transfer rate for evaporation and condensation_______________________________________________________________________

    dmass_EVAP = ZK.mass_dot_evap_funct (HEOS, 
                             c_dict, 
                             n_dict, 
                             L, 
                             A, 
                             g, 
                             pressure_tank, 
                             volume_VAP, 
                             volume_LIQ, 
                             temperature_VAP, 
                             temperature_LIQ, 
                             E = 693,
                             CS = CS_tank)
    
    dmass_COND = ZK.mass_dot_cond_funct (HEOS,
                             pressure_tank, 
                             volume_VAP,
                             temperature_VAP, 
                             temperature_LIQ,
                             E = 693,
                             CS = CS_tank,
                             R_universal = universal_gas_constant,
                             delta_t = time_step)

# Mass transfer rate between the liquid and vapor phases____________________________________________________________________
    dmass_VAP = ZK.mass_dot_VAP_funct (dmass_EVAP, 
                                       dmass_COND, 
                                       dmass_OUT)
    dmass_LIQ = ZK.mass_dot_LIQ_funct (dmass_EVAP, 
                                       dmass_COND)

# Volume change rate for liquid and vapour phases___________________________________________________________________________
    dvolume_VAP = ZK.volume_dot_vap_funct (HEOS, 
                                         pressure_tank, 
                                         mass_VAP, 
                                         mass_LIQ, 
                                         dmass_VAP, 
                                         dmass_LIQ, 
                                         volume_VAP, 
                                         volume_LIQ, 
                                         temperature_VAP, 
                                         temperature_LIQ,
                                         dmass_EVAP, 
                                         dmass_COND,
                                         dmass_OUT,
                                         dheat_in_vap,
                                         dheat_in_liq)
    
    dvolume_LIQ = ZK.volume_dot_liq_funct (dvolume_VAP)
    
# Density change rate for liquid and vapour phases__________________________________________________________________________
    drho_VAP = ZK.density_dot_funct(HEOS, 
                                    mass_VAP, 
                                    dmass_VAP, 
                                    volume_VAP, 
                                    dvolume_VAP)
    
    drho_LIQ = ZK.density_dot_funct(HEOS, 
                                    mass_LIQ, 
                                    dmass_LIQ, 
                                    volume_LIQ, 
                                    dvolume_LIQ)

# Energy transfer rate between the liquid and vapour phases_________________________________________________________________
    denergy_VAP = ZK.energy_dot_vap_funct (HEOS, 
                                           pressure_tank, 
                                           temperature_VAP, 
                                           temperature_LIQ, 
                                           dmass_EVAP, 
                                           dmass_COND, 
                                           dheat_in_vap, 
                                           dvolume_VAP)
    
    denergy_LIQ = ZK.energy_dot_liq_funct (HEOS, 
                                           pressure_tank, 
                                           temperature_LIQ, 
                                           dmass_EVAP, 
                                           dmass_COND, 
                                           dmass_OUT, 
                                           dheat_in_liq, 
                                           dvolume_LIQ, 
                                           threshold = 0.01)
    
# Temperature change rate for the liquid and vapor phases___________________________________________________________________
    dtemperature_VAP = ZK.temperature_dot_funct (HEOS, 
                                                 mass_VAP, 
                                                 dmass_VAP,
                                                 drho_VAP,
                                                 pressure_tank,
                                                 temperature_VAP,
                                                 denergy_VAP)
    
    dtemperature_LIQ = ZK.temperature_dot_funct (HEOS, 
                                                 mass_LIQ, 
                                                 dmass_LIQ, 
                                                 drho_LIQ, 
                                                 pressure_tank, 
                                                 temperature_LIQ, 
                                                 denergy_LIQ)

# Outputting combined system of ODE_________________________________________________________________________________________
    dZdt = np.zeros(12)
    dZdt[0] = dmass_LIQ
    dZdt[1] = dmass_VAP
    dZdt[2] = dmass_COND
    dZdt[3] = dmass_EVAP
    dZdt[4] = dmass_OUT
    dZdt[5] = dtemperature_LIQ
    dZdt[6] = dtemperature_VAP
    dZdt[7] = denergy_LIQ
    dZdt[8] = denergy_VAP
    dZdt[9] = drho_LIQ
    dZdt[10] = drho_VAP
    dZdt[11] = dvolume_LIQ
    dZdt[12] = dvolume_VAP
    
    return dZdt
    
    

# %% Input Data for Testing
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
HEOS.update(CP.PT_INPUTS,56*101325, 298)
print(HEOS.T())
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
V_vap = L_vap*CS_tank
V_liq = L_liq*CS_tank
m_vap = 1
m_liq = 1

# Initial Tank conditions
T_vap = 298.15
T_liq = 298.15
P_tank_initial = 56*ct.one_atm


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


threshold = 0.01
E = 693

# %% Solving System of ODE
ti = 0
tf = 5
step_number = 100
time_step = tf/step_number
t_eval = np.linspace(ti,tf,step_number)
y0 = np.zeros(7)
y0[0] = P_tank_initial
y0[1] = T_vap
y0[2] = T_liq
y0[3] = m_vap
y0[4] = m_liq
y0[5] = V_vap
y0[6] = V_liq


sol = solve_ivp(tank_model_funct,[0,tf],y0,
                method='LSODA',
                t_eval=t_eval, 
                args=(HEOS, threshold, E, CS_tank, c_dict, n_dict, L_tank, CS_tank, g, universal_gas_constant),
                max_step=0.01)
