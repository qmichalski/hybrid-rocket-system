# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 04:41:41 2026

@author: arun2
"""

import CoolProp.CoolProp as CP
import math
import cantera as ct
import numpy as np
import Non_Equilibrium_Tank_Model_ZK as ZK
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# %% Input Data for Testing
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
HEOS_LIQ = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
HEOS_VAP = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
HEOS_S_VAP = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
HEOS_S_LIQ = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
#print(HEOS.T())
# Basic Constants
g = 9.81
universal_gas_constant = 8.314

# Tank Geometry
L_tank = 1.652
R_tank = 19.05/200
CS_tank = math.pi*(R_tank**2)

V_tank = L_tank*CS_tank
#print(V_tank)
# Initial Tank conditions
#T_vap = 298.15
T_liq = 280
P_tank_initial = 56*ct.one_atm

# Initial tank data
fill_percentage = 90/100
L_vap = L_tank*(1 - fill_percentage)
L_liq = L_tank*(fill_percentage)
V_vap = L_vap*CS_tank
V_liq = L_liq*CS_tank

HEOS_VAP.update(CP.PQ_INPUTS,P_tank_initial,1)
rho_vap = HEOS_VAP.rhomass()
m_vap = HEOS_VAP.rhomass()*V_vap
u_vap = HEOS_VAP.umass()

HEOS_LIQ.update(CP.PQ_INPUTS,P_tank_initial,0)
rho_liq = HEOS_LIQ.rhomass()
m_liq = HEOS_LIQ.rhomass()*V_liq
u_liq = HEOS_LIQ.umass()

# Heat Transfer Constants
c_dict = {
    'liquid N2O': 0.15,
    'vapour N2O': 0.15
    }
n_dict = {
    'liquid N2O': 1/3,
    'vapour N2O': 1/3
    }
k_vap = 0.0173
k_liq = 0.15

# Solver Inputs
time_step = 0.01
threshold = 0.01
E = 683

# Injector and outflow data
Cd = 0.55
P_exit = 101325
injector_area = 2.54/1000000

U_vap = u_vap*m_vap
U_liq = u_liq*m_liq
M_tot = m_liq + m_vap
# %% SOLVER INITIALIZATION

ti = 0
tf = 10
step_number = 1000
time_step = tf/step_number
t_eval = np.linspace(ti,tf,step_number)
y0 = np.zeros(6)
y0[0] = m_vap
y0[1] = m_liq
y0[2] = rho_vap
y0[3] = rho_liq
y0[4] = U_vap
y0[5] = U_liq

sol = solve_ivp(ZK.tank_model_funct,[ti,tf],y0,
                method='RK45',
                t_eval=t_eval,
                args=(M_tot, time_step, HEOS, HEOS_VAP, HEOS_LIQ, HEOS_S_VAP, HEOS_S_LIQ, threshold, E, CS_tank, c_dict, n_dict, L_tank, CS_tank, g, k_vap, k_liq, universal_gas_constant, Cd, P_exit, injector_area),
                max_step=time_step,
                rtol=1e99,
                atol=1e99,
                first_step = time_step)
