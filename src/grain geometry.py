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

#all outputs in m, m^2 and m^3
#Functions
#solver for Addapted Finocyl surface area, cross section , volume
def graingeometry_finocyl(geometry_scalar,grainlength,armheight,armwidth,numberofarms,graincentreradius):
    graincentreradiusupdate = graincentreradius + geometry_scalar
    armwidthupdate = (armwidth + (geometry_scalar*2))
    burn_fillet = geometry_scalar
    burn_fillet_area = (burn_fillet**2 - (math.pi*burn_fillet**2))*(numberofarms*2) #accounts for the change in geometery area as the sharp corners burn away
    angle = math.asin((armwidthupdate/2)/graincentreradiusupdate) * 2
    arclength = angle*2*graincentreradiusupdate
    arcarea = ((graincentreradiusupdate**2)*angle/2)-((1/2)*(graincentreradiusupdate**2)*math.sin(angle))
    rectanglesubtraction = (grainlength*arclength*numberofarms)
    boresurfacearea = math.pi*2*graincentreradiusupdate*grainlength
    rectangle_surface_area = (grainlength*(armheight-burn_fillet)*numberofarms*2)+((armwidthupdate-burn_fillet)*grainlength*numberofarms)+(((2 * math.pi* burn_fillet)/4)*grainlength)
    
    if rectanglesubtraction > boresurfacearea and (armheight+graincentreradiusupdate) < chamber_outer_radius :
        
        grainsurfacearea = rectangle_surface_area
        graincrosssection = ((3*(3)**(1/2))/2)*(armwidthupdate**2)+(armheight*armwidthupdate*numberofarms)-burn_fillet_area
        Volumeofgrain = graincrosssection*grainlength
        
    else:
        grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectangle_surface_area
        graincrosssection = math.pi*graincentreradiusupdate**2+(armheight*armwidthupdate*numberofarms)- arcarea*numberofarms - burn_fillet_area
        Volumeofgrain = graincrosssection*grainlength
        
    return(grainsurfacearea,graincrosssection,Volumeofgrain,armheight,armwidthupdate,graincentreradiusupdate)    

        
    """ #dead code may be needed if calculating area after failure point.
    elif(armheight+graincentreradius) > chamber_outer_radius: #grain failure
        armheight = armheight+(chamber_outer_radius-(armheight+graincentreradius))
    """    

#input parameters ignore parameters not used in type of grain
chamber_outer_radius = 51.52/2000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
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
    volumeofgrain = graincrosssection*grainlength
    
    #grain calculation
    final_offset = chamber_outer_radius-graincentreradius
    geometry_scalars = numpy.linspace(0,final_offset,1000)
    
    surface_area_finocyl = []
    grain_cross_section_finocyl = []
    grain_volume_finocyl = []
    tracker = []
    
    for geometry_scalar in geometry_scalars:
        
        grainsurfacearea,graincrosssection,volumeofgrain,armheight,armwidthupdate,graincentreradiusupdate = graingeometry_finocyl(geometry_scalar, grainlength, armheight, armwidth, numberofarms, graincentreradius)
        surface_area_finocyl.append(grainsurfacearea)
        grain_cross_section_finocyl.append(graincrosssection)
        grain_volume_finocyl.append(volumeofgrain)
        tracker.append(armwidthupdate)# for troublshooting
    ''' 
    #could not get to work
    y0 = numpy.zeros(5) # = [0,0]
    y0[0] = numberofarms
    y0[1] = grainlength
    y0[2] = graincentreradius
    y0[3] = armheight
    y0[4] = armwidth
    
    final_offset = chamber_outer_radius-graincentreradius
    geometry_scalar = numpy.linspace(0,final_offset,1000)
    sol = solve_ivp(graingeometry,[0,final_offset],y0,
                    method='RK45',geometry_scalar=geometry_scalar,
                    args=(
                          armheight,
                          armwidth,
                          numberofarms,
                          graincentreradius),
                    max_step=0.001)
        
    print ('grain failure expected at:', final_offset)
    '''
    
    regression_vs_Surface_area_curve = Chebyshev.fit(geometry_scalars ,surface_area_finocyl , deg=6)

    graph.plot(geometry_scalars ,surface_area_finocyl)
    graph.xlabel("regression (m)")
    graph.ylabel("surface area finocyl(m^2)")
    graph.title('regression vs Surface area')
    graph.grid()
    graph.show()
    
    regression_vs_grain_cross_section_curve = Chebyshev.fit(geometry_scalars ,grain_cross_section_finocyl , deg=6)

    graph.plot(geometry_scalars ,grain_cross_section_finocyl)
    graph.xlabel("regression (m)")
    graph.ylabel("grain cross section finocyl(m^2)")
    graph.title('regression vs cross section')
    graph.grid()
    graph.show()
    
    regression_vs_volume_curve = Chebyshev.fit(geometry_scalars ,grain_volume_finocyl , deg=6)

    graph.plot(geometry_scalars ,grain_volume_finocyl)
    graph.xlabel("regression (m)")
    graph.ylabel("grain volume finocyl(m^3)")
    graph.title('regression vs volume') 
    graph.grid()
    graph.show()
    
    print ('grain failure expected at:', final_offset*1000,'mm')
    
               
elif typeofgrain == 'straight bore':
    
    grainsurfacearea = math.pi*2*graincentreradius*grainlength
    graincrosssection = math.pi*2*graincentreradius
    Volumeofgrain = graincrosssection*grainlength

else:
    print ("grain type not implemented yet")


