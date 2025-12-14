# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 16:47:22 2025

@author: quent
"""
from scipy.optimize import root
from scipy.optimize import minimize_scalar
import CoolProp.CoolProp as CP

def _diff_massflux_HEM_choked_specified_Q0_T0(p2,HEOS,h0,s0,fulloutput=False):
    # print(p2)
    HEOS.update(CP.PQ_INPUTS, p2, 0) # saturated liquid
    sL = HEOS.smass()
    hL = HEOS.hmass()
    rhoL = HEOS.rhomass()
    HEOS.update(CP.PQ_INPUTS, p2, 1)
    sV = HEOS.smass()
    hV = HEOS.hmass()
    rhoV = HEOS.rhomass()
    xHEM = (s0-sL)/(sV-sL)
    # print(xHEM)
    hHEM = xHEM*hV+(1-xHEM)*hL
    vHEM = xHEM/rhoV+(1-xHEM)/rhoL
    rhoHEM = 1/vHEM
    GHEM = rhoHEM*(2*(h0-hHEM))**(1/2)
    if fulloutput:
        return(rhoHEM, GHEM/rhoHEM)
    else:
        return(1/GHEM)

def massflux_HEM_choked_specified_Q0_T0(T0, Q0, HEOS):
    HEOS.update(CP.QT_INPUTS, Q0, T0)
    s0 = HEOS.smass()
    h0 = HEOS.hmass()
    p0 = HEOS.p()
    sol = minimize_scalar(_diff_massflux_HEM_choked_specified_Q0_T0,bounds=[p0/2,p0*0.99],args=(HEOS,h0,s0))
    p_choked = sol['x']
    rho_choked,v_choked = _diff_massflux_HEM_choked_specified_Q0_T0(p_choked,HEOS,h0,s0,fulloutput=True)
    return(rho_choked,v_choked,p_choked)

def throat_flow_HEM_Q0_T0(T0, Q0, HEOS, p_downstream):
    rho_choked,v_choked,p_choked = massflux_HEM_choked_specified_Q0_T0(T0, Q0, HEOS)
    HEOS.update(CP.QT_INPUTS, Q0, T0)
    s0 = HEOS.smass()
    h0 = HEOS.hmass()
    p0 = HEOS.p()
    if p_downstream <= p_choked:
        return(rho_choked,v_choked,p_choked)
    else:
        if p_downstream > p0:
            rho0 = HEOS.rhomass()
            return(rho0,0,p0)
        else:
            HEOS.update(CP.PQ_INPUTS, p_downstream, 0)
            sL = HEOS.smass()
            hL = HEOS.hmass()
            rhoL = HEOS.rhomass()
            HEOS.update(CP.PQ_INPUTS, p_downstream, 1)
            sV = HEOS.smass()
            hV = HEOS.hmass()
            rhoV = HEOS.rhomass()
            xHEM = (s0-sL)/(sV-sL)
            # print(xHEM)
            hHEM = xHEM*hV+(1-xHEM)*hL
            vHEM = xHEM/rhoV+(1-xHEM)/rhoL
            rhoHEM = 1/vHEM
            GHEM = rhoHEM*(2*(h0-hHEM))**(1/2)
            return(rhoHEM,GHEM/rhoHEM,p_downstream)

def sound_speed(gas):
    #return gas sound speed dpdrho at constant entropy
    T,P = gas.TP
    rho = gas.density
    P_perturbated = P*(1+1e-5)
    gas.SP = gas.entropy_mass, P_perturbated
    rho_perturbated = gas.density
    gas.TP = T,P
    return(((P-P_perturbated)/(rho-rho_perturbated))**(1/2))

def _diff_choked_flow(P,T0,P0,X0,gas,fulloutput=False):
    gas.TPX = T0,P0,X0
    s0 = gas.entropy_mass
    h0 = gas.enthalpy_mass
    gas.SP = s0,P[0]
    h = gas.enthalpy_mass
    a = sound_speed(gas)
    v = (2*(h0-h))**(1/2)
    if fulloutput:
        return(gas.density,v,P[0])
    else:
        return(a-v)
    
def choked_flow(gas):
    T0, P0, X0 = gas.TPX
    sol = root(_diff_choked_flow,x0=gas.P*0.9,args=(T0,P0,X0,gas))
    P = sol['x']
    rho_choked,v_choked,p_choked = _diff_choked_flow(P,T0,P0,X0,gas,fulloutput=True)
    gas.TPX = T0, P0, X0
    return(rho_choked,v_choked,p_choked)

def throat_flow(gas,p_downstream):
    T0, P0, X0 = gas.TPX
    rho_choked,v_choked,p_choked = choked_flow(gas)
    s0 = gas.entropy_mass
    h0 = gas.enthalpy_mass
    if p_downstream <= p_choked:
        return(rho_choked,v_choked,p_choked)
    else:
        if p_downstream > gas.P:
            return(gas.density,0,gas.P)
        else:
            gas.SP = s0,p_downstream
            rho = gas.density
            h = gas.enthalpy_mass
            v = abs((2*(h0-h))**(1/2))
            gas.TPX=T0, P0, X0
            return(rho,v,p_downstream)

def massflux_DRYER_Q0_T0(T0, Q0, HEOS, p_downstream):
    HEOS.update(CP.QT_INPUTS, Q0, T0)
    s0 = HEOS.smass()
    h0 = HEOS.hmass()
    p0 = HEOS.p()
    rho0 = HEOS.rhomass()
    HEOS.update(CP.PQ_INPUTS, p_downstream, 0)
    sL = HEOS.smass()
    hL = HEOS.hmass()
    rhoL = HEOS.rhomass()
    HEOS.update(CP.PQ_INPUTS, p_downstream, 1)
    sV = HEOS.smass()
    hV = HEOS.hmass()
    rhoV = HEOS.rhomass()
    xHEM = (s0-sL)/(sV-sL)
    # print(xHEM)
    hHEM = xHEM*hV+(1-xHEM)*hL
    vHEM = xHEM/rhoV+(1-xHEM)/rhoL
    rhoHEM = 1/vHEM
    if p0 < p_downstream:
        return(0)
    else:
        GHEM = rhoHEM*(2*(h0-hHEM))**(1/2)
        GSPI = (2*rho0*(p0-p_downstream))**(1/2)
        return((GHEM+GSPI)/2)