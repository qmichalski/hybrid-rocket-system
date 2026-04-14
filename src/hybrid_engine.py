"""
Created on Wed Dec 10 18:05:38 2025

@author: quent,guillaume,arun
"""

import cantera as ct
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import CoolProp.CoolProp as CP
import libCombRegRate
import libMassFlux
import grain_geometry_lib
import Nozzle_Model as NM #please do not abriviate definitions
from mpl_toolkits.mplot3d import Axes3D   # needed for 3D plotting
from scipy.interpolate import griddata

def plot_scattered_surface(x, y, z,
                           grid_size=100,
                           method='linear',
                           show_points=True,
                           cmap='viridis'):


    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    if not (x.shape == y.shape == z.shape):
        raise ValueError("x, y, z must all be 1D with the same length")


    xi = np.linspace(x.min(), x.max(), grid_size)
    yi = np.linspace(y.min(), y.max(), grid_size)
    X, Y = np.meshgrid(xi, yi)   # X, Y are 2D


    Z = griddata((x, y), z, (X, Y), method=method)  # Z is 2D


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, Z, cmap=cmap)

    if show_points:
        ax.scatter(x, y, z, color='red', s=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)

    plt.show()

# MAIN CHAMBER FUNCTION_________________________________________________________________________________________
def combustionDifferentialSystemWithOxidizerTank(t,z,
                                                 section_nozzle,
                                                 ambiant_pressure,
                                                 volume_oxidizer,
                                                 section_injector,
                                                 Y_oxidizer,
                                                 Swr_fun,Pr_fun,combustion_final_radius,grain_length,quench_radius,
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
    # quench_radius = 1e-3 # 1e-3 # mm
    if combustion_radius > combustion_final_radius-quench_radius:
        volume_combustor = Pr_fun(combustion_final_radius-quench_radius)*combustor_length
    else:
        volume_combustor =  Pr_fun (combustion_radius)*combustor_length
    rhoc = mc/volume_combustor # chamber averaged density
    Yc = mYc/mc
    gas.UVY = uc, 1/rhoc, Yc
    pc = gas.P # chamber pressure
    Tc = gas.T # chamber temperature
    gamma = gas.cp/gas.cv
    R = gas.cp - gas.cv
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
    

    # Getting hot gas properties
    mu = gas.viscosity
    print(f"Dynamic Viscosity: {mu:.4e} Pa*s")
    
    # Pr = mu * Cp / lambda
    cp = gas.cp_mass  # Specific heat at constant pressure
    lambda_val = gas.thermal_conductivity
    Pr = mu * cp / lambda_val
    print(f"Prandtl Number: {cp} ")

    print('t:{:4.2f}s|pc:{:4.2f} bar|Tc:{:4.0f}K|mdoto:{:1.3f}|r/rmax:{:03.0f}'.format(t,pc/1e5,Tc,mdot_oxidizer,combustion_radius/combustion_final_radius*100))
    if combustion_radius < combustion_final_radius-quench_radius and mdot_oxidizer > 0: # grain run out controller
        area_combustion = Swr_fun(combustion_radius)
        section_combustor = Pr_fun(combustion_radius)
        mdot_fuel,T_burnt_products,Y_burnt_products,h_reactants_mix,regression_rate = libCombRegRate.solveCombustion(
                          h_fuel,Y_fuel,pc,h_oxidizer,Y_oxidizer,
                          T_abs,cp_abs,rho_abs,hv,area_combustion,
                          mdot_oxidizer,section_combustor,combustor_length,combustion_eff,gas)
    else:
        area_combustion = 0
        section_combustor = Pr_fun(combustion_final_radius-quench_radius)
        mdot_fuel = 0
        #T_burnt_products = To
        Y_burnt_products = Y_oxidizer
        h_reactants_mix = h_oxidizer
        regression_rate = 0
    combustionVolumedVdt = area_combustion*regression_rate # to calculate expansion work
    mdot_reactants = mdot_oxidizer+mdot_fuel
    dmcdt = mdot_reactants-massflow_throat
    dUdct = mdot_reactants*h_reactants_mix-massflow_throat*hc-pc*combustionVolumedVdt
    dmYcdt = mdot_reactants*Y_burnt_products - massflow_throat*Yc
    dmodt = -mdot_oxidizer
    dUdot = -mdot_oxidizer*ho
    
    dzdt = np.zeros(len(z))
    
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
               po,To,xo,rhoo,mdot_oxidizer,mdot_fuel,gamma,R,Pr,mu)
    else:
        return(dzdt)

# COMBUSTION QUENCHING__________________________________________________________________________________________
# def combustion_quenching(t,z,section_nozzle,
#                         ambiant_pressure,
#                         volume_oxidizer,
#                         section_injector,
#                         Y_oxidizer,
#                         Swr_fun,Pr_fun,combustion_final_radius,grain_length,quench_radius,
#                         Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
#                         gas,HEOS):
#     return(combustion_final_radius-z[0])

def tank_empty(t,z,section_nozzle,
                        ambiant_pressure,
                        volume_oxidizer,
                        section_injector,
                        Y_oxidizer,
                        Swr_fun,Pr_fun,combustion_final_radius,grain_length,quench_radius,
                        Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                        gas,HEOS):
    combustion_radius = z[0] # combustion radius
    mc = z[1] # combustor mass
    Uc = z[2] # combustor internal energy
    mo = z[3] # mass of oxidizer in tank
    Uo = z[4] # oxidizer internal energy
    mYc = z[5:]
    uo = Uo/mo # chamber averaged specific internal energy
    rhoo = mo/volume_oxidizer # chamber averaged density
    HEOS.update(CP.DmassUmass_INPUTS, rhoo, uo)
    po = HEOS.p()
    uc = Uc/mc # chamber averaged specific internal energy
    combustor_length = grain_length
    # quench_radius = 1e-3 # mm
    if combustion_radius > combustion_final_radius-quench_radius:
        volume_combustor = Pr_fun(combustion_final_radius-quench_radius)*combustor_length
    else:
        volume_combustor =  Pr_fun (combustion_radius)*combustor_length
    rhoc = mc/volume_combustor # chamber averaged density
    Yc = mYc/mc
    gas.UVY = uc, 1/rhoc, Yc
    pc = gas.P # chamber pressure
    return(po-pc-10) #this is the tollerence

# combustion_quenching.terminal = True
# combustion_quenching.direction = -1
tank_empty.terminal = True
tank_empty.direction = -1

# speciesToKeep = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'CO', 'CO2', 'NO', 'N2', 'N2O', 'C3H8']
# speciesToKeep = ['O2', 'N2', 'N2O', 'C3H8']
# species = [gas.species(name) for name in speciesToKeep]
# create the new reduced mechanism
# gas = ct.Solution(thermo='ideal-gas',species=species)

plt.rcParams['axes.grid'] = True
#fluid libraries
HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
gas = ct.Solution('gri30_highT.yaml')

#input parameters for grain ignore parameters not used in type of grain in m ___________________________________

chamber_outer_radius = 25.76/1000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
typeofgrain = 'stargrain'
numberofarms = 8 #only used for grains with radial features
grain_length = 358/1000
graincentreradius = 10/1000
armheight = 8/1000 #only used for grains with radial features
armwidth = 4.229/1000 #only used for grains with radial features

# FUEL DATA_____________________________________________________________________________________________________

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

ambiant_pressure = 101235
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
#diameter_nozzle = 2*9.43e-3 # m
#section_nozzle = np.pi*(diameter_nozzle/2)**2
volume_oxidizer = 4.63e-3 # in m^3
section_injector = 18.85e-6 #np.pi*(10e-3/2)**2
rho_abs = 1080 # kg/m3
cp_abs = 1500 # J/kg/K
hv = 1.8e6 # J/kg
combustion_eff = 0.8
quench_radius = 1e-3

Swr_fun,Pr_fun,v_fun,combustion_final_radius,grain_length = grain_geometry_lib.grain_solver(chamber_outer_radius,typeofgrain,numberofarms,grain_length,graincentreradius,armheight,armwidth)


# %% Solver Setup

throat_radius_initial = 8
throat_radius_final = 10
throat_radius_spacing = 1
exit_radius_initial = 13
exit_radius_final = 18
Exit_radius_spacing = 1
R_T = range(throat_radius_initial, throat_radius_final, throat_radius_spacing+1)
R_E = range(exit_radius_initial, exit_radius_final, Exit_radius_spacing+1)

T_curve_area = []
R_T_choice = []
R_E_choice = []

for throat_radius in R_T:
    # NOZZLE DATA___________________________________________________________________________________________________
    
    
    #throat_radius = 9.43 # In mm
    #exit_radius = 17.5 # In mm
    throat_area = math.pi * throat_radius ** 2 / 10**6
    Area_throat = throat_area
    section_nozzle = Area_throat
    
    """
    M_design = 2.79
    Area_exit = math.pi * 19.00 ** 2 / 10**6
    throat_area = math.pi * 8.89 ** 2 / 10**6
    # %%
    """
    
    # SOLUTION______________________________________________________________________________________________________
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
    tf = 8 #this cannot be less then the run time. Set it to 20 it will self terminate at the oxidiser run out condition.
    t_eval = np.linspace(0,tf,1000)
    
    sol = solve_ivp(combustionDifferentialSystemWithOxidizerTank,[0,tf],y0,
                    method='LSODA',t_eval=t_eval,events=(tank_empty), 
                    args=(
                    section_nozzle,
                    ambiant_pressure,
                    volume_oxidizer,
                    section_injector,
                    Y_oxidizer,
                    Swr_fun,Pr_fun,combustion_final_radius,grain_length,quench_radius,
                    Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                    gas,HEOS),
                    max_step=0.01)
    
#%% exit radius variation loop
    for exit_radius in R_E:
        # %%EXTRACTING DATA_______________________________________________________________________________________________
        pcs = []
        Tcs = []
        pos = []
        Tos = []
        xos = []
        mdoto = []
        mdot_total = []
        mdotf = []
        thrust = []
        TVA = [] # time varience authority
        MU = []
        PR = []
        mdottotf = []
        C_STAR = []
        pcs_G = []
        stop = 0
        thrust_sum = 0
        mass_flow_rate_sum = 0
        gamma = 1.21
        
        for i,t in enumerate(sol['t']):
            # try:
            z = np.zeros(len(y0))
            z[0] = sol['y'][0,i]
            z[1] = sol['y'][1,i]
            z[2] = sol['y'][2,i]
            z[3] = sol['y'][3,i]
            z[4:] = sol['y'][4:,i] 
            pc,Tc,uc,rhoc,Yc,rho_throat, v_throat, p_throat,po,To,xo,rhoo,mdot_oxidizer,mdot_fuel,gamma,R, Pr_val, mu_val = combustionDifferentialSystemWithOxidizerTank(t,z,
                                                                                                                     section_nozzle,
                                                                                                                     ambiant_pressure,
                                                                                                                     volume_oxidizer,
                                                                                                                     section_injector,
                                                                                                                     Y_oxidizer,
                                                                                                                     Swr_fun,Pr_fun,combustion_final_radius,grain_length,quench_radius,
                                                                                                                     Y_fuel,T_abs,cp_abs,rho_abs,hv,combustion_eff,
                                                                                                                     gas,HEOS,fulloutput=True)
            pcs.append(pc/1e3)
            pcs_G.append(pc**((1-gamma)/(2*gamma)))
            Tcs.append(Tc)
            pos.append(po/1e3)
            Tos.append(To)
            xos.append(xo)
            MU.append(mu_val)
            PR.append(Pr_val)
            mdoto.append(mdot_oxidizer*1e3)
            mdotf.append(mdot_fuel*1e3)
            
            
            c_star = Area_throat*pc/(mdot_oxidizer+mdot_fuel)
            C_STAR.append(c_star)

            
            if mdot_fuel > 0:
                mdottotf.append(mdot_oxidizer/mdot_fuel)
            else:
                mdottotf.append(0)
            try:
                thrust_nozzle, P_throat, T_throat, P_exit, T_exit, v_exit, a_exit, m_dot_exit = NM.nozzle_sol_funct(1*ct.one_atm, 100, pc, Tc, "CO2:17.03, H2O:9.45, N2:0.5", throat_radius/1000, exit_radius/1000)
            except:
                thrust_nozzle, P_throat, T_throat, P_exit, T_exit, v_exit, a_exit, m_dot_exit = 0,0,0,0,0,0,0,0
            thrust.append(thrust_nozzle)
            #T_exit = Tc/ (1 + (gamma-1)*(M_design**2)/2)
            #a_exit = math.sqrt(gamma*R*T_exit)
            #V_exit = M_design*a_exit
            #m_dot_exit = rho_throat*throat_area*v_throat
            #P_exit = pc /( (1 + (gamma-1)*(M_design**2)/2)**(gamma/(gamma-1)))
            if thrust_nozzle > 0:
                thrust[i] = thrust_nozzle
            else:
                if stop < 2:
                    thrust_finish = i
                    thrust_time = t
                    stop =stop + 1
                thrust[i] = 0
            thrust_sum = thrust_sum +thrust[i]
            mass_flow_rate_sum = mass_flow_rate_sum + m_dot_exit
            TVA.append(t)
            # except:
                # print("oh no", t)
                # tf = t
                # # t_eval = np.linspace(0,tf,100)
                # break
            mdot_total.append(m_dot_exit)
        fuel_mass = rho_abs*(Pr_fun(combustion_final_radius)-Pr_fun(sol['y'][0]))*grain_length
        
        
        # %%PLOTTING______________________________________________________________________________________________________
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
        # sol['t'] = TVA
        ax[0].plot(TVA,Tcs,label='T combustion chamber',color='red')
        ax[0].plot(TVA,pcs,label='P combustion chamber',color='black')
        ax[0].plot(TVA,pos,label='P N2O tank',color='blue')
        ax[1].plot(TVA,mdoto,label='$\dot{m}_{N2O}$',color='blue')
        ax[1].plot(TVA,mdotf,label='$\dot{m}_{ABS}$',color='black')
        ax[2].plot(TVA,mdottotf,color='black')
        ax[3].plot(TVA,fuel_mass*1e3,color='black')
        ax[0].legend()
        ax[1].legend()
        ax[-1].set_xlabel('Time,s')
        ax[0].set_ylabel('T, K | P, kPa')
        ax[1].set_ylabel('$\dot{m}, g/s$')
        ax[2].set_ylabel('Mixture Ratio')
        ax[3].set_ylabel('Grain mass, g')
        plt.show()
        
        max_thrust = max(thrust)
        xmax = thrust.index(max_thrust)
        if thrust_finish ==0:
            thrust_finish = 1
        average_thrust = sum(thrust)/thrust_finish
        max_thrust_text = ('peak thrust '+ str(round(max_thrust))+ ' N')
        average_thrust_text = ('average thrust '+ str(round(average_thrust))+ ' N')
        average_thrust_display = average_thrust+10
        display_point = (max(TVA)/1.5)
        total_impulse = average_thrust*thrust_time
        
        plt.annotate(max_thrust_text, xy = (xmax,max_thrust))
        plt.annotate(average_thrust_text, xy = (display_point,average_thrust_display))
        plt.hlines(average_thrust, 0, max(TVA),colors = 'blue')
        plt.ylim([0,max(thrust)+100])
        plt.plot(TVA,thrust,color = 'red')
        plt.xlabel('Time, s')
        plt.ylabel('Thrust, N')
        plt.show()
        
        
        plt.plot(TVA, C_STAR,color = 'Blue')
        plt.xlabel('Time, s')
        plt.ylabel('Characteristic Velocity, m/s')
        plt.show()
        
        plt.plot(TVA, mdot_total,color = 'Blue')
        plt.xlabel('Time, s')
        plt.ylabel('Mass flow rate, Kg/s')
        plt.show()
        
        #print(total_impulse)
        impulse = np.trapezoid(thrust, x=TVA)
        #print('Area under P_c curve = ',np.trapezoid(pcs, x=TVA)*1000)
        #print('Area under P_c^G curve = ',np.trapezoid(pcs_G , x=TVA))
        print('Area under Thrust curve = ',impulse)
        R_T_choice.append(throat_radius)
        R_E_choice.append(exit_radius)
        T_curve_area.append(impulse)
        print(throat_radius,exit_radius)
"""
fig = plt.figure()
ax = plt.axes(projection = '3d')
surf = ax.plot_surface(R_T_choice, R_E_choice, T_curve_area)
ax.set_title('Impulse Surface Plot')
ax.set_xlabel('Throat Radius (mm)')
ax.set_ylabel('Exit Radius (mm)')
ax.set_zlabel('Total Impulse (Ns)')
"""
#plot_scattered_surface(R_T_choice,R_E_choice,T_curve_area)
R_T_choice = np.array(R_T_choice)
R_E_choice = np.array(R_E_choice)
T_curve_area = np.array(T_curve_area)
max_value = T_curve_area.max()
max_index = T_curve_area.argmax()

print("Best throat radius = ",R_T_choice[max_index])
print("Best exit radius   = ",R_E_choice[max_index])
print("Best Impulse       = ",max_value)
