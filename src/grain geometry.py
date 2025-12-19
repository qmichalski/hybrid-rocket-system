# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 15:04:45 2025

@author: Guillaume Knight
"""

import cantera as ct
import math
import matplotlib.pyplot as graph
import numpy
import time
from numpy.polynomial import Chebyshev
import CoolProp.CoolProp as coolprop
from CoolProp.CoolProp import PropsSI as propsi
from scipy.integrate import solve_ivp
mech = 'Mevel2017.yaml'

#Functions stored at end
#solver for Addapted Finocyl
def graingeometry(geometry_scalar,grainlength,armheight,armwidth,numberofarms,graincentreradius):
    graincentreradius = graincentreradius+geometry_scalar
    armheight = armheight + geometry_scalar
    armwidth = armwidth + geometry_scalar
    burn_fillet = geometry_scalar
    burn_fillet_area = (burn_fillet**2 - (math.pi*burn_fillet**2))*(numberofarms*2) #accounts for the change in geometery area as the sharp corners burn away
    angle = math.asin((armwidth/2)/graincentreradius) * 2
    arclength = angle*2*graincentreradius
    arcarea = ((graincentreradius**2)*angle/2)-((1/2)*(graincentreradius**2)*math.sin(angle))
    rectanglesubtraction = grainlength*arclength*numberofarms
    boresurfacearea = math.pi*2*graincentreradius*grainlength
    rectanglearea = (grainlength*(armheight-burn_fillet)*numberofarms*2)+((armwidth-burn_fillet)*grainlength*numberofarms)+(((2 * math.pi* burn_fillet)/4)*grainlength)
    
    if rectanglesubtraction > boresurfacearea and (armheight+graincentreradius) < chamber_outer_radius :
        
        grainsurfacearea = rectanglearea
        graincrosssection = ((3*(3)**(1/2))/2)*(armwidth**2)+(armheight*armwidth*numberofarms)-burn_fillet_area
        Volumeofgrain = graincrosssection*grainlength
        
    elif(armwidth*grainlength*numberofarms) > chamber_outer_radius:
        
    else:
        grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectanglearea
        graincrosssection = math.pi*graincentreradius**2+(armheight*armwidth*numberofarms)- arcarea*numberofarms - burn_fillet_area
        Volumeofgrain = graincrosssection*grainlength
    return(grainsurfacearea,graincrosssection,Volumeofgrain,armheight,armwidth,graincentreradius)    

#input parameters ignore parameters not used in type of grain
chamber_outer_radius = 51.52/1000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
typeofgrain = 'Addapted Finocyl'
numberofarms = 6 #only used for grains with radial features
grainlength = 358/1000
graincentreradius = 10/1000
armheight = 8/1000 #only used for grains with radial features
armwidth = 4.229/1000 #only used for grains with radial features

#start conditions
if typeofgrain == 'Addapted Finocyl':
    
    angle = math.asin((armwidth/2)/graincentreradius) * 2
    arclength = angle*2*graincentreradius
    arcarea = ((graincentreradius**2)*angle/2)-((1/2)*(graincentreradius**2)*math.sin(angle))
    rectanglesubtraction = grainlength*arclength*numberofarms
    boresurfacearea = math.pi*2*graincentreradius*grainlength
    rectanglearea = (grainlength*armheight*numberofarms*2)+(armwidth*grainlength*numberofarms)
    grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectanglearea
    graincrosssection = math.pi*graincentreradius**2+(armheight*armwidth*numberofarms)- arcarea*6
    Volumeofgrain = graincrosssection*grainlength       
               
elif typeofgrain == 'straight bore':
    
    grainsurfacearea = math.pi*2*graincentreradius*grainlength
    graincrosssection = math.pi*2*graincentreradius
    Volumeofgrain = graincrosssection*grainlength

else:
    print ("grain type not implemented yet")
    
final_offset = chamber_outer_diamter-graincentreradius
geometry_scalar = numpy.linspace(0,final_offset,1000)
sol = solve_ivp(graingeometry,[0,final_offset]
                method='LSODA',geometry_scalar=geometry_scalar,
                args=(chamber_pressure,
                volume_oxidizer,
                section_injector,
                HEOS),
                max_step=0.001)
    
'''

statetemp = 0
y0 = numpy.zeros(2) # = [0,0]
masstank = HEOS.rhomass()*capacitytank
internalenergybymass = masstank*HEOS.umass()
y0[0] = masstank
y0[1] = internalenergybymass
totaltime = 10 # s
timestep = numpy.linspace(0,totaltime,1000)
sol = solve_ivp(massflowtotal,[0,totaltime],y0,
                method='LSODA',timestep=timestep,
                args=(temperaturen2o,
                        capacitytank,
                        areainjector,
                        HEOS,
                        graincrosssection,
                        armheight,
                        armwidth,
                        graincentreradius,
                        ),
                        max_step=0.001)
pos = numpy.zeros(len(timestep))
Tos = numpy.zeros(len(timestep))
xos = numpy.zeros(len(timestep))
mdoto = numpy.zeros(len(timestep))
mos = sol['y'][0,:]
for i,t in enumerate(timestep):
    z = numpy.zeros(len(y0))
    z[0] = sol['y'][0,i]
    z[1] = sol['y'][1,i]
    # z[2:] = sol['y'][2:,i] 
    pressuren2o ,temperaturen2o ,fractionvapour ,averagedensity ,oxidisermassflow = massflowtotal(time,z,
                                                        temperaturen2o,
                                                        pressurechamber,
                                                        capacitytank,
                                                        areainjector,
                                                        HEOS,
                                                        graincrosssection,
                                                        statetemp,
                                                        armheight,
                                                        armwidth,
                                                        graincentreradius)
    pos[i] = pressuren2o
    Tos[i] = temperaturen2o
    xos[i] = fractionvapour
    mdoto[i] = oxidisermassflow
    
graph.plot(sol['t'],Tos-273.15,label='T')
graph.plot(sol['t'],pos/1e5,label='P')
graph.xlabel('time,s')
graph.ylabel('T, K | P, kPa')
graph.show()

graph.plot(sol['t'],mdoto)
graph.xlabel('time,s')
graph.ylabel('massflow oxidizer, kg/s')
graph.show()

graph.plot(sol['t'][1:],xos[1:]*mos[1:],label='Gaseous mass')
graph.plot(sol['t'][1:],(1-xos[1:])*mos[1:],label='Liquid mass')
graph.legend()
graph.xlabel('time,s')
graph.ylabel('mass, kg')
graph.show()

graph.plot(pressuren2o,)
graph.ylabel("mass (kg)")
graph.xlabel("time (s)")
graph.grid()
graph.title('tank discharge')
graph.show()

#def-diff-massflux(P-)

#dmcdt = masstank
#dUdct = *h_reactants_mix-oxidisermassflow*hc
#dmYcdt = mdot_reactants*Y_burnt_products - massflow_throat*Yc
'''
