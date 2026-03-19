import CoolProp.CoolProp as CP
import math
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
# %% Test Function
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

HEOS_VAP.update(CP.PT_INPUTS,P_tank_initial,298)
rho_vap = HEOS_VAP.rhomass()
m_vap = HEOS_VAP.rhomass()*V_vap
u_vap = HEOS_VAP.umass()

HEOS_LIQ.update(CP.PT_INPUTS,P_tank_initial,280)
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
E = 693
A = CS_tank
# Injector and outflow data
Cd = 0.55
P_exit = 101325
injector_area = 2.54/1000000

U_vap = u_vap*m_vap
U_liq = u_liq*m_liq
M_tot = m_liq + m_vap

# %% Functions
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')

P = 56*ct.one_atm
T = 298
m = 1.8

L_VAP_SURF = V_vap/(2*A)
L_LIQ_SURF = L_tank/2 - L_VAP_SURF
HEOS_S_LIQ.update(CP.PQ_INPUTS,P_tank_initial,0)
HEOS_S_VAP.update(CP.PQ_INPUTS,P_tank_initial,1)
temperature_SURF = HEOS_S_LIQ.T()
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

dmass_COND = mass_dot_cond_funct (HEOS,
                     HEOS_VAP,
                     P_tank_initial,
                     V_vap,
                     universal_gas_constant,
                     time_step
                     )

print(dheat_liq_to_surf)
print(dheat_surf_to_vap)
print(dmass_COND)

for T in range(270,298):
    HEOS.update(CP.QT_INPUTS, 1, T)
    P_sat = HEOS.p()
    plt.plot(T,P_sat/ct.one_atm,'.')
plt.grid()


"""
    print(HEOS.hmass(),'J/Kg')
    print(HEOS.rhomass(),'Kg/m^3')
    print('U = ',HEOS.umass()*m,'J')
    print('U_calc = ', m*HEOS.cpmass()*T, 'J')
"""