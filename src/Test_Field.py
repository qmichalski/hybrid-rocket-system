import CoolProp.CoolProp as CP
import math
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# %% Input Data for Testing
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')

P = 56*ct.one_atm
T = 298
m = 1.8

for P in range(1,60):
    HEOS.update(CP.PT_INPUTS, P*ct.one_atm, 280)
    rho = HEOS.rhomass()
    plt.plot(P,rho,'.')
plt.grid()
"""
    print(HEOS.hmass(),'J/Kg')
    print(HEOS.rhomass(),'Kg/m^3')
    print('U = ',HEOS.umass()*m,'J')
    print('U_calc = ', m*HEOS.cpmass()*T, 'J')
"""