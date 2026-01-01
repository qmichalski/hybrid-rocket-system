# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 21:56:13 2025

@author: quent
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import CoolProp.CoolProp as CP
import math

def regRateEilers(T_combustion,T_abs,Pr,cp_abs,rho_abs,hv,mdot_ox,viscosity_ox,section_combustor,combustor_length):
    # Shannon D. Eilers∗ and Stephen A. Whitmore, 2008, DOI: 10.2514/1.33804
    # Correlation of Hybrid Rocket Propellant Regression Measurements with Enthalpy-Balance Model Predictions
    r = 0.047/(Pr**(2/3)*rho_abs)*(cp_abs*(T_combustion-T_abs)/hv)**(0.23)
    r = r*(mdot_ox/section_combustor)**(4/5)
    r = r*(viscosity_ox/combustor_length)**(1/5)
    return(r) # m/s

def _diff_regRate(mdot_fuel,h_fuel,Y_fuel,pc,h_oxidizer,Y_oxidizer,
                  T_abs,cp_abs,rho_abs,hv,area_combustion,
                  mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas,fulloutput=False,choice='Eilers'):
    mdot_reactants = mdot_oxidizer+mdot_fuel
    gas.HPY = h_oxidizer, pc, Y_oxidizer
    viscosity_ox = gas.viscosity # Pa.s
    h_reactants_mix = (mdot_oxidizer*h_oxidizer+mdot_fuel*h_fuel)/(mdot_reactants) 
    Y_reactants_mix = (mdot_oxidizer*Y_oxidizer+mdot_fuel*Y_fuel)/(mdot_reactants)
    gas.HPY = h_reactants_mix, pc, Y_reactants_mix 
    gas.equilibrate('HP')
    Y_adiabatic_burnt_products = gas.Y
    Y_burnt_products = combustion_eff*Y_adiabatic_burnt_products+(1-combustion_eff)*Y_reactants_mix
    gas.HPY = h_reactants_mix, pc, Y_burnt_products 
    if choice == 'Eilers':
        Pr = gas.viscosity*gas.cp_mass/gas.thermal_conductivity
        T_burnt_products = gas.T
        r = regRateEilers(T_burnt_products,T_abs,Pr,cp_abs,rho_abs,hv,mdot_oxidizer,viscosity_ox,section_combustor,combustor_length)
    mdot_fuel_calc = rho_abs*r*area_combustion
    # print(mdot_fuel[0],mdot_fuel_calc,r,mdot_oxidizer/mdot_fuel[0])
    if fulloutput:
        return(T_burnt_products,Y_burnt_products,h_reactants_mix,r)
    return(mdot_fuel-mdot_fuel_calc)

def solveCombustion(h_fuel,Y_fuel,Pc,h_oxidizer,Y_oxidizer,
                  T_abs,cp_abs,rho_abs,hv,area_combustion,
                  mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas):
    # getting stoichiometric mixture ratio for further calculation
    gas.set_equivalence_ratio(1,fuel=Y_fuel,oxidizer=Y_oxidizer,basis='mass')
    MRst = sum(gas.Y[np.where(Y_oxidizer>0)[0]])/sum(gas.Y[np.where(Y_fuel>0)[0]])
    sol = root_scalar(_diff_regRate,x0=mdot_oxidizer/MRst,bracket=(0,mdot_oxidizer*1e3),
               args=(h_fuel,Y_fuel,Pc,h_oxidizer,Y_oxidizer,
                     T_abs,cp_abs,rho_abs,hv,area_combustion,
                     mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas),method="brentq")
    mdot_fuel = sol.root# ['x']
    T_burnt_products,Y_burnt_products,h_reactants_mix,r = _diff_regRate(
                     mdot_fuel,h_fuel,Y_fuel,Pc,h_oxidizer,Y_oxidizer,
                     T_abs,cp_abs,rho_abs,hv,area_combustion,
                     mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas,fulloutput=True)
    return(mdot_fuel,T_burnt_products,Y_burnt_products,h_reactants_mix,r)
"""
def combustionDifferentialSystem(t,z,Swr_fun,Pr_fun,combustion_final_radius,grain_length,
                                 T_abs,cp_abs,rho_abs,hv,
                                 section_nozzle,
                                 ambiant_pressure,
                                 volume_oxidizer,
                                 section_injector,
                                 gas,fulloutput=False):
    combustion_radius = z
    
    mdot_oxidizer = 0.3
    # mdot_fuel = mdot_oxidizer/9 # calculate from regression rate
    # rho_throat, v_throat, p_throat = throat_flow(gas,ambiant_pressure)
    # massflow_throat = rho_throat*v_throat*section_nozzle
    #calculating burnt gas properties
    gas.TPX = 300,1e5,'N2O:1'
    h_oxidizer = gas.enthalpy_mass
    Y_oxidizer = gas.Y
    gas.TPX = 300,1e5,'ABS:1'
    h_fuel = gas.enthalpy_mass
    Y_fuel = gas.Y 
    Pc = 20e5
    area_combustion = Swr_fun(combustion_radius)
    section_combustor = Pr_fun(combustion_radius)
    combustor_length = grain_length
    mdot_fuel,T_burnt_products,Y_burnt_products,h_reactants_mix,regression_rate = solveCombustion(
                      h_fuel,Y_fuel,Pc,h_oxidizer,Y_oxidizer,
                      T_abs,cp_abs,rho_abs,hv,area_combustion,
                      mdot_oxidizer,section_combustor,combustor_length,gas)
    
    dzdt = regression_rate
    if fulloutput:
        return(mdot_fuel,regression_rate)
    else:
        return(dzdt)
    
if __name__ == "__main__":

    def combustion_quenching(t,z,Swr_fun,Pr_fun,combustion_final_radius,grain_length,
                                     T_abs,cp_abs,rho_abs,hv,
                                     section_nozzle,
                                     ambiant_pressure,
                                     volume_oxidizer,
                                     section_injector,
                                     gas,HEOS):
        return(combustion_final_radius-z[0])

    combustion_quenching.terminal = True
    combustion_quenching.direction = -1    

    gas = ct.Solution('gri30_highT.yaml')
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
    # iFuel = gasfull.species_index('ABS')
    # pviscosityFuel = gasfull.get_viscosity_polynomial(gasfull.species_index('CO2'))
    # gasfull.set_viscosity_polynomial(iFuel,pviscosityFuel)
    # speciesToKeep = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'CO', 'CO2', 'NO', 'N2', 'N2O', 'C3H8', 'ABS']
    # gas = ct.Solution(thermo='ideal-gas',species=species)
    # gas = gasfull['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'CO', 'CO2', 'NO', 'N2', 'N2O', 'C3H8', 'ABS']
    # gas.TPX = 300,5e5,'O2:0.21,N2:0.79'
    T0 = 300
    T_abs = T0
    
    gas.TPY = T0,1e5,'ABS:1'
    Y_fuel = gas.Y
    h_fuel = gas.enthalpy_mass
    gas.TPY = T0,1e5,'N2O:1'
    Y_oxidizer = gas.Y
    h_oxidizer = gas.enthalpy_mass
    Pc = 25e5
    combustor_diameter = 15e-3
    combustor_length = 250e-3 
    area_combustion = np.pi*combustor_diameter*combustor_length
    section_combustor = np.pi*(combustor_diameter/2)**2
    volume_combustor = section_combustor*combustor_length#784058e-9#
    diameter_nozzle = 2*9.43e-3 # m
    section_nozzle = np.pi*(diameter_nozzle/2)**2
    ambiant_pressure = 1e5
    volume_oxidizer = 4.63e-3
    section_injector = 18.85e-6 #np.pi*(10e-3/2)**2
    rho_abs = 1080 # kg/m3
    cp_abs = 1500 # J/kg/K
    hv = 1.8e6 # J/kg
    mdot_oxidizer = 0.4 # kg/s
    combustion_eff = 0.8
    
    mdot_fuel,T_burnt_products,Y_burnt_products,h_reactants_mix,r = solveCombustion(
                      h_fuel,Y_fuel,Pc,h_oxidizer,Y_oxidizer,
                      T_abs,cp_abs,rho_abs,hv,area_combustion,
                      mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas)
    
    print(T_burnt_products)
    # Swr_fun,Pr_fun,combustion_final_radius,grain_length = combustionBurningGrainGeometry(grainType='Circular_1')
    
    # tf = 16
    # t_eval = np.linspace(0,tf,1000)
    # y0 = [0]
    # sol = solve_ivp(combustionDifferentialSystem,[0,tf],y0,
    #                 method='LSODA',t_eval=t_eval,events=combustion_quenching,
    #                 args=(Swr_fun,Pr_fun,combustion_final_radius,grain_length,
    #                                                  T_abs,cp_abs,rho_abs,hv,
    #                                                  section_nozzle,
    #                                                  ambiant_pressure,
    #                                                  volume_oxidizer,
    #                                                  section_injector,
    #                                                  gas),
    #                 max_step=0.1)
    
    # plt.plot(sol['t'],sol['y'][0],'.')
    # print(sol['t'][-1],sol['y'][0][-1]/combustion_final_radius)
    # mdot_fuels = np.linspace(mdot_oxidizer/20,mdot_oxidizer,20)
    # diffs = np.zeros(len(mdot_fuels))
    # for i,mdot_fuel in enumerate(mdot_fuels):
    #     diff = _diff_regRate([mdot_fuel],h_fuel,Y_fuel,Pc,h_oxidizer,Y_oxidizer,
    #                       T_abs,cp_abs,rho_abs,hv,area_combustion,
    #                       mdot_oxidizer,section_combustor,combustor_length,gas,fulloutput=False)
    #     diffs[i] = diff
    # plt.plot(mdot_oxidizer/mdot_fuels,diffs)
    """