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
#solver for Addapted Finocyl surface area, cross section , volume
def graingeometry(geometry_scalar,grainlength,armheight,armwidth,numberofarms,graincentreradius):
    graincentreradius = graincentreradius + geometry_scalar
    armwidth = armwidth + geometry_scalar
    burn_fillet = geometry_scalar
    burn_fillet_area = (burn_fillet**2 - (math.pi*burn_fillet**2))*(numberofarms*2) #accounts for the change in geometery area as the sharp corners burn away
    angle = math.asin((armwidth/2)/graincentreradius) * 2
    arclength = angle*2*graincentreradius
    arcarea = ((graincentreradius**2)*angle/2)-((1/2)*(graincentreradius**2)*math.sin(angle))
    rectanglesubtraction = grainlength*arclength*numberofarms
    boresurfacearea = math.pi*2*graincentreradius*grainlength
    rectangle_surface_area = (grainlength*(armheight-burn_fillet)*numberofarms*2)+((armwidth-burn_fillet)*grainlength*numberofarms)+(((2 * math.pi* burn_fillet)/4)*grainlength)
    
    if rectanglesubtraction > boresurfacearea and (armheight+graincentreradius) < chamber_outer_radius :
        
        grainsurfacearea = rectangle_surface_area
        graincrosssection = ((3*(3)**(1/2))/2)*(armwidth**2)+(armheight*armwidth*numberofarms)-burn_fillet_area
        Volumeofgrain = graincrosssection*grainlength
    else:
        grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectangle_surface_area
        graincrosssection = math.pi*graincentreradius**2+(armheight*armwidth*numberofarms)- arcarea*numberofarms - burn_fillet_area
        Volumeofgrain = graincrosssection*grainlength
        
    return(grainsurfacearea,graincrosssection,Volumeofgrain,armheight,armwidth,graincentreradius)    

        
    """ #dead code may be needed if calculating area after failure point.
    elif(armheight+graincentreradius) > chamber_outer_radius: #grain failure
        armheight = armheight+(chamber_outer_radius-(armheight+graincentreradius))
    """    

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
    rectangle_surface_area = (grainlength*armheight*numberofarms*2)+(armwidth*grainlength*numberofarms)
    grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectangle_surface_area
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
sol = solve_ivp(graingeometry,[0,final_offset],
                method='LSODA',geometry_scalar=geometry_scalar,
                args=(chamber_pressure,
                volume_oxidizer,
                section_injector,
                HEOS),
                max_step=0.001)
    

