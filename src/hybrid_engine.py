"""
Created on Wed Dec 10 18:05:38 2025

@author: quent,guillaume
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import CoolProp.CoolProp as CP
import libCombRegRate
import libMassFlux
import grain_geometry_lib

#fluid libraries
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
gas = ct.Solution('gri30_highT.yaml')

#input parameters for grain ignore parameters not used in type of grain in m

chamber_outer_radius = 25.76/1000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
typeofgrain = 'Addapted Finocyl'
numberofarms = 6 #only used for grains with radial features
grain_length = 358/1000
graincentreradius = 10/1000
armheight = 8/1000 #only used for grains with radial features
armwidth = 4.229/1000 #only used for grains with radial features

#inital fluid conditions

fuel_name = 'ABS'
fuel_composition = 'C:17.03, H:18.9, N:1'
fuel_eof = 90.312*1e6# J/mol
fuel = ct.Species(fuel_name,fuel_composition)
fuel.thermo = ct.ConstantCp(200,5000,101325,
                              (298.15,
                              fuel_eof,# J/kmol
                              0,
                              0))
tran = ct.GasTransportData()
tran.set_customary_units('nonlinear', 3.75, 141.40, 0.0, 2.60, 13.00)
fuel.transport = tran
gas.add_species(fuel)

ambiant_pressure = 1e5
T0 = 300
oxidizer = 'N2O:1'
gas.TPX = T0,ambiant_pressure,oxidizer
Y_oxidizer = gas.Y
fuel = 'ABS:1'
gas.TPX = T0,ambiant_pressure,fuel
Y_fuel = gas.Y
gas.TPX = T0,ambiant_pressure,'O2:0.21,N2:0.79'
# gas.TPY = T0,1e5,'N2O:3,C3H8:1'
HEOS.update(CP.QT_INPUTS, 0, T0)
# combustor_length = 150e-3 
T_abs = T0
diameter_nozzle = 2*9.43e-3 # m
section_nozzle = np.pi*(diameter_nozzle/2)**2
volume_oxidizer = 4.63e-3
section_injector = 18.85e-6 #np.pi*(10e-3/2)**2
rho_abs = 1080 # kg/m3
cp_abs = 1500 # J/kg/K
hv = 1.8e6 # J/kg
combustion_eff = 0.8

#program start

Swr_fun,Pr_fun,v_fun,combustion_final_radius,grain_length = grain_geometry_lib.grain_solver(chamber_outer_radius,typeofgrain,numberofarms,grain_length,graincentreradius,armheight,armwidth)

def combustionDifferentialSystemWithOxidizerTank(t,z,
                                                 section_nozzle,
                                                 ambiant_pressure,
                                                 volume_oxidizer,
                                                 section_injector,
                                                 Y_oxidizer,
                                                 Swr_fun,Pr_fun,combustion_final_radius,grain_length,
                                                 Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                                                 gas,HEOS,fulloutput=False):
    combustion_radius = z[0] # combustion radius
    mc = z[1] # combustor mass
    Uc = z[2] # combustor internal energy
    mo = z[3] # mass of oxidizer in tank
    Uo = z[4] # oxidizer internal energy
    mYc = z[5:]
    uo = Uo/mo # chamber averaged specific internal energy
    rhoo = mo/volume_oxidizer # chamber averaged density
    HEOS.update(CP.DmassUmass_INPUTS, rhoo, uo)
    To = HEOS.T()
    ho = HEOS.hmass()
    uc = Uc/mc # chamber averaged specific internal energy
    combustor_length = grain_length
    volume_combustor =  Pr_fun(combustion_radius)*combustor_length
    rhoc = mc/volume_combustor # chamber averaged density
    Yc = mYc/mc
    gas.UVY = uc, 1/rhoc, Yc
    pc = gas.P # chamber pressure
    Tc = gas.T # chamber temperature
    hc = gas.enthalpy_mass # chamber average enthalpy
    # if t >= 0.1: # and t <= 1.1:
    massflux_injector = libMassFlux.massflux_DRYER_Q0_T0(To, 0, HEOS, pc)
    mdot_oxidizer = 0.61*massflux_injector*section_injector # Cd to be passed as coefficient
    rho_throat, v_throat, p_throat = libMassFlux.throat_flow(gas,ambiant_pressure)
    massflow_throat = rho_throat*v_throat*section_nozzle
    #calculating burnt gas properties
    gas.TPY = To,1e5,Y_oxidizer
    h_oxidizer = gas.enthalpy_mass
    gas.TPY = 300,1e5,Y_fuel
    h_fuel = gas.enthalpy_mass
    area_combustion = Swr_fun(combustion_radius)
    section_combustor = Pr_fun(combustion_radius)
    mdot_fuel,T_burnt_products,Y_burnt_products,h_reactants_mix,regression_rate = libCombRegRate.solveCombustion(
                      h_fuel,Y_fuel,pc,h_oxidizer,Y_oxidizer,
                      T_abs,cp_abs,rho_abs,hv,area_combustion,
                      mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas)
    mdot_reactants = mdot_oxidizer+mdot_fuel
    dmcdt = mdot_reactants-massflow_throat
    dUdct = mdot_reactants*h_reactants_mix-massflow_throat*hc
    dmYcdt = mdot_reactants*Y_burnt_products - massflow_throat*Yc
    dmodt = -mdot_oxidizer
    dUdot = -mdot_oxidizer*ho
    
    dzdt = np.zeros(len(z))
    print('t:{:4.2f}s|pc:{:4.2f} bar|Tc:{:4.0f}K|mdott:{:1.3f}|Tp:{:03.0f}K'.format(t,pc/1e5,Tc,massflow_throat,T_burnt_products))
    # print('O2:{:4.2f}s|pc:{:4.2f} bar|Tc:{:4.0f}K|mdott:{:1.3f}kgs'.format(t,pc/1e5,Tc,massflow_throat))
    dzdt[0] = regression_rate
    dzdt[1] = dmcdt
    dzdt[2] = dUdct
    dzdt[3] = dmodt
    dzdt[4] = dUdot
    dzdt[5:] = dmYcdt
    if fulloutput:
        HEOS.update(CP.DmassUmass_INPUTS, rhoo, uo)
        po = HEOS.p()
        xo = HEOS.Q()
        return(pc,Tc,uc,rhoc,Yc,rho_throat, v_throat, p_throat,
               po,To,xo,rhoo,mdot_oxidizer,mdot_fuel)
    else:
        return(dzdt)

def combustion_quenching(t,z,section_nozzle,
                        ambiant_pressure,
                        volume_oxidizer,
                        section_injector,
                        Y_oxidizer,
                        Swr_fun,Pr_fun,combustion_final_radius,grain_length,
                        Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                        gas,HEOS):
    return(combustion_final_radius-z[0])

combustion_quenching.terminal = True
combustion_quenching.direction = -1

# speciesToKeep = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'CO', 'CO2', 'NO', 'N2', 'N2O', 'C3H8']
# speciesToKeep = ['O2', 'N2', 'N2O', 'C3H8']
# species = [gas.species(name) for name in speciesToKeep]
# create the new reduced mechanism
# gas = ct.Solution(thermo='ideal-gas',species=species)


mo0 = HEOS.rhomass()*volume_oxidizer
Uo0 = mo0*HEOS.umass()
volume_combustor = Pr_fun(0)*grain_length
mc0 = gas.density*volume_combustor
Uc0 = mc0*gas.int_energy_mass
mYc0 = mc0*gas.Y
y0 = np.zeros((5+len(mYc0)))
y0[0] = 0
y0[1] = mc0
y0[2] = Uc0
y0[3] = mo0
y0[4] = Uo0
y0[5:] = mYc0
tf = 10
t_eval = np.linspace(0,tf,100)
sol = solve_ivp(combustionDifferentialSystemWithOxidizerTank,[0,tf],y0,
                method='LSODA',t_eval=t_eval,events=combustion_quenching, 
                args=(
                section_nozzle,
                ambiant_pressure,
                volume_oxidizer,
                section_injector,
                Y_oxidizer,
                Swr_fun,Pr_fun,combustion_final_radius,grain_length,
                Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                gas,HEOS),
                max_step=0.1)

pcs = np.zeros(len(t_eval))
Tcs = np.zeros(len(t_eval))
pos = np.zeros(len(t_eval))
Tos = np.zeros(len(t_eval))
xos = np.zeros(len(t_eval))
mdoto = np.zeros(len(t_eval))
mdotf = np.zeros(len(t_eval))
for i,t in enumerate(t_eval):
    z = np.zeros(len(y0))
    z[0] = sol['y'][0,i]
    z[1] = sol['y'][1,i]
    z[2] = sol['y'][2,i]
    z[3] = sol['y'][3,i]
    z[4:] = sol['y'][4:,i] 
    pc,Tc,uc,rhoc,Yc,rho_throat,v_throat,p_throat,po,To,xo,rhoo,mdot_oxidizer,mdot_fuel = combustionDifferentialSystemWithOxidizerTank(t,z,
                                                                                                             section_nozzle,
                                                                                                             ambiant_pressure,
                                                                                                             volume_oxidizer,
                                                                                                             section_injector,
                                                                                                             Y_oxidizer,
                                                                                                             Swr_fun,Pr_fun,combustion_final_radius,grain_length,
                                                                                                             Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                                                                                                             gas,HEOS,fulloutput=True)
    pcs[i] = pc
    Tcs[i] = Tc
    pos[i] = po
    Tos[i] = To
    xos[i] = xo
    mdoto[i] = mdot_oxidizer
    mdotf[i] = mdot_fuel

fuel_mass = rho_abs*(Pr_fun(combustion_final_radius)-Pr_fun(sol['y'][0]))*grain_length

inchToMM = 25.4
combustionAndFlameMax = 88 # mm
plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams["font.size"] = 12
plt.rcParams['mathtext.fontset'] = "cm"
plt.rcParams['mathtext.rm'] = "cm"
plt.rcParams['mathtext.it'] = "cm"
plt.rcParams['mathtext.bf'] = "cm"
fig1 = plt.figure(figsize=(combustionAndFlameMax/inchToMM*2, combustionAndFlameMax/inchToMM*2),dpi=400)
gs1 = fig1.add_gridspec(4, 1, hspace=0.1,height_ratios=[3,1,1,1])
ax = gs1.subplots(sharex=True)
ax[0].plot(sol['t'],Tcs,label='T combustion chamber',color='red')
ax[0].plot(sol['t'],pcs/1e3,label='P combustion chamber',color='black')
ax[0].plot(sol['t'],pos/1e3,label='P N2O tank',color='blue')
ax[1].plot(sol['t'],mdoto*1e3,label='$\dot{m}_{N2O}$',color='blue')
ax[1].plot(sol['t'],mdotf*1e3,label='$\dot{m}_{ABS}$',color='black')
ax[2].plot(sol['t'],mdoto/mdotf,color='black')
ax[3].plot(sol['t'],fuel_mass*1e3,color='black')
ax[0].legend()
ax[1].legend()
ax[-1].set_xlabel('Time,s')
ax[0].set_ylabel('T, K | P, kPa')
ax[1].set_ylabel('$\dot{m}, g/s$')
ax[2].set_ylabel('Mixture Ratio')
ax[3].set_ylabel('Grain mass, g')
# ax[-1].set_xlim([])
# ax[-1].set_xscale('log')