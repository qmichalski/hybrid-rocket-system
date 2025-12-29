# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 15:04:45 2025

@author: Guillaume Knight
"""

import math
import matplotlib.pyplot as graph
import numpy
from numpy.polynomial import Chebyshev
from scipy.interpolate import interp1d
import sys

'''
chamber_outer_radius = 25.76/1000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
typeofgrain = 'Straight Bore'
numberofarms = 6 #only used for grains with radial features
grain_length = 358/1000
graincentreradius = 10/1000
armheight = 8/1000 #only used for grains with radial features
armwidth = 4.229/1000 #only used for grains with radial features
'''

#all outputs in m, m^2 and m^3
#solver for Addapted Finocyl surface area, cross section , volume

#Functions=

def grain_geometry_finocyl(geometry_scalar,grain_length,armheight,armwidth,numberofarms,graincentreradius,chamber_outer_radius):
    
    graincentreradiusupdate = graincentreradius + geometry_scalar
    a=(((graincentreradius+geometry_scalar)**2)-(armwidth/2)**2)**(1/2)-((graincentreradius**2-(armwidth/2)**2))**(1/2)
    b=(((graincentreradius+geometry_scalar)**2)-((armwidth/2)**2))**(1/2)-(((graincentreradius+geometry_scalar)**2) - ((armwidth+(geometry_scalar*2))/2)**2)**(1/2)
    c = a-b
    armheightupdate = (armheight+geometry_scalar)- c
    armwidthupdate = (armwidth + (geometry_scalar*2))
    burn_fillet = geometry_scalar
    burn_fillet_area = ((((burn_fillet*2)**2) - (math.pi*burn_fillet**2))/4)*(numberofarms*2) #accounts for the change in geometery area as the sharp corners burn away
    angle = math.asin((armwidthupdate/2)/graincentreradiusupdate) * 2
    arclength = angle*graincentreradiusupdate
    arcarea = ((graincentreradiusupdate**2)*angle/2)-((1/2)*(graincentreradiusupdate**2)*math.sin(angle))
    rectanglesubtraction = (grain_length*arclength*numberofarms)
    boresurfacearea = math.pi*2*graincentreradiusupdate*grain_length
    rectangle_surface_area = (grain_length*(armheightupdate-burn_fillet)*numberofarms*2)+((armwidthupdate-(burn_fillet*2))*grain_length*numberofarms)+((((2 * math.pi* burn_fillet)/4)*grain_length)*numberofarms*2)

    if  (boresurfacearea - rectanglesubtraction) <  0 and (armheight+graincentreradiusupdate) < chamber_outer_radius :
        offset = (((graincentreradius)**2)-(armwidth/2)**2)**(1/2)
        armheightupdate = armheight-((armwidthupdate/2)/(math.tan(0.5235987756))-offset)+geometry_scalar
        rectangle_surface_area = (grain_length*(armheightupdate-burn_fillet)*numberofarms*2)+((armwidthupdate-(burn_fillet*2))*grain_length*numberofarms)+((((2 * math.pi* burn_fillet)/4)*grain_length)*numberofarms*2)
        grainsurfacearea = rectangle_surface_area
        graincrosssection = ((3*((3)**(1/2)))/2)*(armwidthupdate**2)+(armheightupdate*armwidthupdate*numberofarms)-(burn_fillet_area)
        volumeofgrain = graincrosssection*grain_length
        
    else:
        grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectangle_surface_area
        graincrosssection = math.pi*graincentreradiusupdate**2+(armheightupdate*armwidthupdate*numberofarms) - (arcarea*numberofarms + burn_fillet_area)
        volumeofgrain = graincrosssection*grain_length
        
    return(grainsurfacearea,graincrosssection,volumeofgrain,armheight,armwidthupdate,graincentreradiusupdate,rectanglesubtraction,)    

        
    """ #dead code may be needed if calculating area after failure point.
    elif(armheight+graincentreradius) > chamber_outer_radius: #grain failure
        armheight = armheight+(chamber_outer_radius-(armheight+graincentreradius))
    """
#start conditions
def grain_solver(chamber_outer_radius,typeofgrain,numberofarms,grain_length,graincentreradius,armheight,armwidth):
    if numberofarms == 0:
        
        typeofgrain = 'Straight Bore'
        
    if typeofgrain == 'Addapted Finocyl':
        
            angle = math.asin((armwidth/2)/graincentreradius) * 2
            arclength = angle*graincentreradius
            arcarea = ((graincentreradius**2)*angle/2)-((1/2)*(graincentreradius**2)*math.sin(angle))
            rectanglesubtraction = grain_length*arclength*numberofarms
            boresurfacearea = math.pi*2*graincentreradius*grain_length
            rectangle_surface_area = (grain_length*armheight*numberofarms*2)+(armwidth*grain_length*numberofarms)
            grainsurfacearea = (boresurfacearea-rectanglesubtraction)+rectangle_surface_area
            graincrosssection = math.pi*graincentreradius**2+(armheight*armwidth*numberofarms)- arcarea*numberofarms
            volumeofgrain = graincrosssection*grain_length
            #print (rectanglesubtraction)
            #print (boresurfacearea-rectanglesubtraction)
            #print (grainsurfacearea)
            
            #visulise inputed grain
            x = []
            y = []
            x1 = []
            y1 = []
            
            x,y,x1,y1 = graingeometry_finocyl_plotter_inital_shape(grain_length,armheight,armwidth,numberofarms,graincentreradius,chamber_outer_radius)
            graph.plot(x,y)
            graph.plot((x),(y))
            graph.xlabel("x (mm)")
            graph.ylabel("y (mm)")
            graph.title('grain inital geometry')
            graph.grid()
            graph.show()
            x.clear()
            y.clear()
            #grain calculation
            final_offset = chamber_outer_radius-(graincentreradius+armheight)
            combustion_final_radius = final_offset
            geometry_scalars = numpy.linspace(0,final_offset,50) #adjusts resolution of plot
            
            surface_area_finocyl = []
            grain_cross_section_finocyl = []
            grain_volume_finocyl = []
            tracker = []
            
    
            for geometry_scalar in geometry_scalars:
    
                x,y = grain_geometry_finocyl_plotter(geometry_scalar,grain_length,armheight,armwidth,numberofarms,graincentreradius,chamber_outer_radius)
                    
                graph.plot(x,y)
            
            graph.xlabel("x (mm)")
            graph.ylabel("y (mm)")
            graph.title('regression profile')
            graph.grid()
            graph.show()
            
            geometry_scalars = numpy.linspace(0,final_offset,1000)
            for geometry_scalar in geometry_scalars:
                
                grainsurfacearea,graincrosssection,volumeofgrain,armheight,armwidthupdate,graincentreradiusupdate,rectanglesubtraction = grain_geometry_finocyl (geometry_scalar, grain_length, armheight, armwidth, numberofarms, graincentreradius, chamber_outer_radius)
                surface_area_finocyl.append(grainsurfacearea)
                grain_cross_section_finocyl.append(graincrosssection)
                grain_volume_finocyl.append(volumeofgrain)
                tracker.append(rectanglesubtraction)# for troublshooting
            
            regression_vs_Surface_area_curve = Chebyshev.fit(geometry_scalars ,surface_area_finocyl , deg=6)
            Swr_fun = regression_vs_Surface_area_curve
            
            graph.plot(geometry_scalars ,surface_area_finocyl)
            graph.xlabel("regression (m)")
            graph.ylabel("surface area finocyl(m^2)")
            graph.title('regression vs surface area')
            graph.show()
            
            regression_vs_grain_cross_section_curve = Chebyshev.fit(geometry_scalars ,grain_cross_section_finocyl , deg=6)
            Pr_fun = regression_vs_grain_cross_section_curve
        
            graph.plot(geometry_scalars ,grain_cross_section_finocyl)
            graph.xlabel("regression (m)")
            graph.ylabel("grain cross section finocyl(m^2)")
            graph.title('regression vs cross section')
            graph.show()
            
            regression_vs_volume_curve = Chebyshev.fit(geometry_scalars ,grain_volume_finocyl , deg=6)
            v_fun = regression_vs_volume_curve
        
            graph.plot(geometry_scalars ,grain_volume_finocyl)
            graph.xlabel("regression (m)")
            graph.ylabel("grain volume finocyl(m^3)")
            graph.title('regression vs volume') 
            graph.show()
            
            print ('grain failure expected at:', final_offset*1000,'mm of regressed material')
                
    elif typeofgrain == 'Straight Bore':
        N = 1000
        Swr_discretes = numpy.linspace(2*math.pi*graincentreradius,2*math.pi*chamber_outer_radius,N)*grain_length
        Swr_discretes = numpy.concatenate([Swr_discretes,[0]])
        Pr_discretes = numpy.linspace(math.pi*graincentreradius**2,math.pi*chamber_outer_radius**2,N)
        Pr_discretes = numpy.concatenate([Pr_discretes,[Pr_discretes[-1]]])
        v_discretes = numpy.linspace(math.pi*graincentreradius**2,math.pi*chamber_outer_radius**2,N)*grain_length
        v_discretes = numpy.concatenate([v_discretes,[v_discretes[-1]]])
        radii_discretes = numpy.linspace(0,chamber_outer_radius-graincentreradius,N)
        radii_discretes = numpy.concatenate([radii_discretes,[radii_discretes[-1]]])
        Swr_fun = interp1d(radii_discretes,Swr_discretes) # wet surface function
        Pr_fun = interp1d(radii_discretes,Pr_discretes) # port section
        v_fun = interp1d(radii_discretes,v_discretes)
        combustion_final_radius = chamber_outer_radius-graincentreradius
        graph.plot(radii_discretes ,Swr_discretes)
        graph.xlabel("regression (m)")
        graph.ylabel("surface area")
        graph.title('regression vs surface area circle bore (m^2)')
        graph.show()
    
        graph.plot(radii_discretes ,Pr_discretes)
        graph.xlabel("regression (m)")
        graph.ylabel("grain cross section")
        graph.title('regression vs cross section circle bore (m^2)')
        graph.show()
    
        graph.plot(radii_discretes ,v_discretes)
        graph.xlabel("regression (m)")
        graph.ylabel("grain volume")
        graph.title('regression vs volume circle bore(m^3)')
        graph.show()
        
    
            
    
    else:
        print ("grain type not implemented yet")
    return(Swr_fun,Pr_fun,v_fun,combustion_final_radius,grain_length)
#grain_solver(chamber_outer_radius,typeofgrain,numberofarms,grain_length,graincentreradius,armheight,armwidth)

def grain_geometry_finocyl_plotter (geometry_scalar,grainlength,armheight,armwidth,numberofarms,graincentreradius,chamber_outer_radius):
   #plot start
   geometry_scalar = geometry_scalar*1000
   chamber_outer_radius = chamber_outer_radius*1000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
   graincentreradius = graincentreradius*1000
   armheight = armheight*1000 #only used for grains with radial features
   armwidth = armwidth*1000 #only used for grains with radial features
   graincentreradiusupdate = graincentreradius + geometry_scalar
   a=(((graincentreradius+geometry_scalar)**2)-(armwidth/2)**2)**(1/2)-((graincentreradius**2-(armwidth/2)**2))**(1/2)
   b=(((graincentreradius+geometry_scalar)**2)-((armwidth/2)**2))**(1/2)-(((graincentreradius+geometry_scalar)**2) - ((armwidth+(geometry_scalar*2))/2)**2)**(1/2)
   c = a-b
   armheightupdate = (armheight+geometry_scalar)- c
   armwidthupdate = (armwidth + (geometry_scalar*2))
   burn_fillet = geometry_scalar
   angle = math.asin((armwidthupdate/2)/graincentreradiusupdate) * 2
   anglebore = ((2*math.pi)-(angle*numberofarms))/numberofarms
   arclength = angle*graincentreradiusupdate
   circumference_bore = graincentreradiusupdate*2*math.pi
   circumference_arc = arclength*numberofarms
   angle_between_arms = (2*math.pi)/numberofarms

   x = []
   y = []
   x1 = []
   y1 = []

   radiuscurve = []
   theacurve = []
   radiuscurve2 =  []
   theacurve2 = []

   pitwo = 2*math.pi
   outercurve = numpy.linspace(0, pitwo ,500)
   for thea in outercurve:
       x1current = math.sin(thea)*chamber_outer_radius
       y1current = math.cos(thea)*chamber_outer_radius
       x1.append(x1current), y1.append(y1current)
   graph.plot((x1) ,(y1))
   graph.gca().set_aspect(1.0)
   ax = graph.gca()
   ax.set_xlim([-chamber_outer_radius, chamber_outer_radius])
   ax.set_ylim([-chamber_outer_radius, chamber_outer_radius])
   if (armheight+graincentreradiusupdate) > chamber_outer_radius:
       finalregression = geometry_scalar
       return(finalregression)
       sys.exit()
   if  circumference_bore - circumference_arc <  0 and (armheight+graincentreradiusupdate) < chamber_outer_radius:
       armtrack = 0
       run =  0
       while armtrack < numberofarms:
           offset = (((graincentreradius)**2)-(armwidth/2)**2)**(1/2)
           armheightupdate = armheight-((armwidthupdate/2)/(math.tan(0.5235987756))-offset)+geometry_scalar
           armtrack = armtrack + 1
           
           thea = (armtrack - 1) * angle_between_arms
           offsetstart = graincentreradius -((graincentreradius**2)-((armwidth/2)**2))**(1/2)
           startheight = armheight+graincentreradius-offsetstart
           startheight = startheight + geometry_scalar
           xstart = math.cos(thea)*(startheight)
           ystart = math.sin(thea)*(startheight)
           x.append(xstart), y.append(ystart)
          
           thea = ((armtrack - 1) * angle_between_arms) + math.pi/2
           startlength  = ((armwidthupdate)/2)-burn_fillet
           yrun = math.sin(thea)*(startlength)
           xrise = math.cos(thea)*(startlength)
           xendpoint = xstart+xrise
           yendpoint = ystart+yrun
           startarmx = numpy.linspace(xstart,xendpoint,100)
           startarmy = numpy.linspace(ystart,yendpoint,100) #breaks the arm into 100 points
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
           xcurrent = xpoint
           ycurrent = ypoint
           ypoint2 =ypoint - armwidthupdate+burn_fillet*2
           
           #burn radius calculations
           if run == 0:
               curve = numpy.linspace(0, burn_fillet ,500)
               for curvey in curve:
                   curvex = ((burn_fillet**2)-(curvey**2))**(1/2)
                   ystore = ypoint+curvey
                   xstore = startheight-(burn_fillet-curvex)
                   r = ((xstore**2)+(ystore**2))**(1/2)
                   theac = math.acos(xstore/r)
                   radiuscurve.append(r)
                   theacurve.append(theac)
                   
               curve2 = numpy.linspace(burn_fillet,0 ,500)
               for curvey in curve2:
                   curvex = ((burn_fillet**2)-(curvey**2))**(1/2)
                   ystore = ypoint2-curvey
                   xstore = startheight-(burn_fillet-curvex)
                   r = ((xstore**2)+(ystore**2))**(1/2)
                   theac = math.asin(ystore/r)
                   radiuscurve2.append(r)
                   theacurve2.append(theac)
                   run = 1
                   
           curvey = 0
           #burn fillet
           for i in range (0, 500):
               angle_update = (armtrack - 1) * angle_between_arms
               thea2 = theacurve[i] + angle_update
               xcurrent = math.cos(thea2)*radiuscurve[i]
               ycurrent = math.sin(thea2)*radiuscurve[i]
               x.append(xcurrent)
               y.append(ycurrent)
           #drop for arm        
           thea1 = (angle/2) 
           thea = (armtrack - 1) * angle_between_arms
           drop = armheightupdate-burn_fillet
           yrun = math.sin(thea)*(drop)
           xrise = math.cos(thea)*(drop)
           xendpoint = xcurrent-xrise
           yendpoint = ycurrent-yrun
           startarmx = numpy.linspace(xcurrent,xendpoint,100)
           startarmy = numpy.linspace(ycurrent,yendpoint,100)
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
           xcurrent = xpoint
           ycurrent = ypoint
          #next arm begins
           thea = ((armtrack) * angle_between_arms) 
           rise = (armheightupdate-burn_fillet)
           yrun = math.sin(thea)*(rise)
           xrise = math.cos(thea)*(rise)
           xendpoint = xcurrent+xrise
           yendpoint = ycurrent+yrun
           startarmx = numpy.linspace(xcurrent,xendpoint,100)
           startarmy = numpy.linspace(ycurrent,yendpoint,100)
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
           xcurrent = xpoint
           ycurrent = ypoint  
           
           for i in range (0, 500):
                   angle_update = (armtrack) * angle_between_arms
                   thea2 = theacurve2[i] + angle_update
                   xcurrent = math.cos(thea2)*radiuscurve2[i]
                   ycurrent = math.sin(thea2)*radiuscurve2[i]
                   x.append(xcurrent)
                   y.append(ycurrent)
           
           thea1 = (angle/2)
           thea = ((armtrack) * angle_between_arms)+ math.pi/2
           startlength  = ((armwidthupdate)/2)-burn_fillet
           yrun = math.sin(thea)*(startlength)
           xrise = math.cos(thea)*(startlength)
           xendpoint = xcurrent+xrise
           yendpoint = ycurrent+yrun
           startarmx = numpy.linspace(xcurrent,xendpoint,100)
           startarmy = numpy.linspace(ycurrent,yendpoint,100) #breaks the arm into 100 points
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
       return x,y
   else:
       armtrack = 0
       run = 0
       #sin and cos can be inverted to rotate grain
       while armtrack < numberofarms:
           armtrack = armtrack + 1
           thea = (armtrack - 1) * angle_between_arms
           offsetstart = graincentreradiusupdate-((graincentreradiusupdate**2)-((armwidthupdate/2)**2))**(1/2)
           startheight = armheightupdate+graincentreradiusupdate-offsetstart
           xstart = math.cos(thea)*(startheight)
           ystart = math.sin(thea)*(startheight)
           x.append(xstart), y.append(ystart)
          
           thea = ((armtrack - 1) * angle_between_arms) + math.pi/2
           startlength  = ((armwidthupdate)/2)-burn_fillet
           yrun = math.sin(thea)*(startlength)
           xrise = math.cos(thea)*(startlength)
           xendpoint = xstart+xrise
           yendpoint = ystart+yrun
           startarmx = numpy.linspace(xstart,xendpoint,100)
           startarmy = numpy.linspace(ystart,yendpoint,100) #breaks the arm into 100 points
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
           xcurrent = xpoint
           ycurrent = ypoint
           ypoint2 =ypoint - armwidthupdate+burn_fillet*2
           
           #burn radius calculations
           if run == 0:
               curve = numpy.linspace(0, burn_fillet ,500)
               for curvey in curve:
                   curvex = ((burn_fillet**2)-(curvey**2))**(1/2)
                   ystore = ypoint+curvey
                   xstore = startheight-(burn_fillet-curvex)
                   r = ((xstore**2)+(ystore**2))**(1/2)
                   theac = math.acos(xstore/r)
                   radiuscurve.append(r)
                   theacurve.append(theac)
                   
               curve2 = numpy.linspace(burn_fillet,0 ,500)
               for curvey in curve2:
                   curvex = ((burn_fillet**2)-(curvey**2))**(1/2)
                   ystore = ypoint2-curvey
                   xstore = startheight-(burn_fillet-curvex)
                   r = ((xstore**2)+(ystore**2))**(1/2)
                   theac = math.asin(ystore/r)
                   radiuscurve2.append(r)
                   theacurve2.append(theac)
                   run = 1
                   
           curvey = 0
           #burn fillet when needed
           if burn_fillet > 0:
               for i in range (0, 500):
                   angle_update = (armtrack - 1) * angle_between_arms
                   thea2 = theacurve[i] + angle_update
                   xcurrent = math.cos(thea2)*radiuscurve[i]
                   ycurrent = math.sin(thea2)*radiuscurve[i]
                   x.append(xcurrent)
                   y.append(ycurrent)
           #drop for arm        
           thea1 = (angle/2) 
           thea = (armtrack - 1) * angle_between_arms
           drop = armheightupdate-burn_fillet
           yrun = math.sin(thea)*(drop)
           xrise = math.cos(thea)*(drop)
           xendpoint = xcurrent-xrise
           yendpoint = ycurrent-yrun
           startarmx = numpy.linspace(xcurrent,xendpoint,100)
           startarmy = numpy.linspace(ycurrent,yendpoint,100)
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
           xcurrent = xpoint
           ycurrent = ypoint
           
           #curve of the bore
           thea1 =(angle/2)+(angle_between_arms*(armtrack-1))
           thea2 = thea1 + anglebore  #calculates the angle of the curve
           curvethea = numpy.linspace(thea1,thea2,500)
           for thea in curvethea:
               xcurrent = math.cos(thea)*graincentreradiusupdate
               ycurrent = math.sin(thea)*graincentreradiusupdate
               x.append(xcurrent), y.append(ycurrent)
           #angled second arm
           
           thea = ((armtrack) * angle_between_arms) 
           rise = (armheightupdate-burn_fillet)
           yrun = math.sin(thea)*(rise)
           xrise = math.cos(thea)*(rise)
           xendpoint = xcurrent+xrise
           yendpoint = ycurrent+yrun
           startarmx = numpy.linspace(xcurrent,xendpoint,100)
           startarmy = numpy.linspace(ycurrent,yendpoint,100)
           for xpoint in startarmx:
               x.append(xpoint) 
       
           for ypoint in startarmy:
               y.append(ypoint)
           xcurrent = xpoint
           ycurrent = ypoint  
           if burn_fillet > 0:
               for i in range (0, 500):
                   angle_update = (armtrack) * angle_between_arms
                   thea2 = theacurve2[i] + angle_update
                   xcurrent = math.cos(thea2)*radiuscurve2[i]
                   ycurrent = math.sin(thea2)*radiuscurve2[i]
                   x.append(xcurrent)
                   y.append(ycurrent)
           
           thea1 = (angle/2)
           thea = ((armtrack) * angle_between_arms)+ math.pi/2
           startlength  = ((armwidthupdate)/2)-burn_fillet
           yrun = math.sin(thea)*(startlength)
           xrise = math.cos(thea)*(startlength)
           xendpoint = xcurrent+xrise
           yendpoint = ycurrent+yrun
           startarmx = numpy.linspace(xcurrent,xendpoint,100)
           startarmy = numpy.linspace(ycurrent,yendpoint,100) #breaks the arm into 100 points
           for xpoint in startarmx:
               x.append(xpoint) 
           for ypoint in startarmy:
               y.append(ypoint)
       return x,y     
   

def graingeometry_finocyl_plotter_inital_shape (grainlength,armheight,armwidth,numberofarms,graincentreradius,chamber_outer_radius,):
   import math
   import numpy
   import matplotlib.pyplot as graph
   #plot start
   chamber_outer_radius = chamber_outer_radius*1000 #unless for some strange reason you are making a pressure vessel out of a non round cross section.
   graincentreradius = graincentreradius*1000
   armheight = armheight*1000 #only used for grains with radial features
   armwidth = armwidth*1000 #only used for grains with radial features
   geometry_scalar = 0
   graincentreradiusupdate = graincentreradius + geometry_scalar
   a=(((graincentreradius+geometry_scalar)**2)-(armwidth/2)**2)**(1/2)-((graincentreradius**2-(armwidth/2)**2))**(1/2)
   b=(((graincentreradius+geometry_scalar)**2)-((armwidth/2)**2))**(1/2)-(((graincentreradius+geometry_scalar)**2) - ((armwidth+(geometry_scalar*2))/2)**2)**(1/2)
   c = a-b
   armheightupdate = (armheight+geometry_scalar)- c
   armwidthupdate = (armwidth + (geometry_scalar*2))
   burn_fillet = geometry_scalar
   angle = math.asin((armwidthupdate/2)/graincentreradiusupdate) * 2
   arclength = angle*graincentreradiusupdate
   anglebore = ((2*math.pi)-(angle*numberofarms))/numberofarms
   angle_between_arms = (2*math.pi)/numberofarms
   circumference_bore = graincentreradiusupdate*2*math.pi
   circumference_arc = arclength*numberofarms

   x = []
   y = []
   x1 = []
   y1 = []
   radiuscurve = []
   theacurve = []
   radiuscurve2 =  []
   theacurve2 = []

   pitwo = 2*math.pi
   outercurve = numpy.linspace(0, pitwo ,500)
   for thea in outercurve:
       x1current = math.sin(thea)*chamber_outer_radius
       y1current = math.cos(thea)*chamber_outer_radius
       x1.append(x1current), y1.append(y1current)
   graph.plot((x1) ,(y1))
   graph.gca().set_aspect(1.0)
   ax = graph.gca()
   ax.set_xlim([-chamber_outer_radius, chamber_outer_radius])
   ax.set_ylim([-chamber_outer_radius, chamber_outer_radius])
   armtrack = 0
   run = 0
   if  circumference_bore - circumference_arc <  0 and geometry_scalar == 0 :
       print("shape not possible")
       print("I don't get paid enough for this crap")
       quit()
   else:
       while armtrack < numberofarms:
               armtrack = armtrack + 1
               thea = (armtrack - 1) * angle_between_arms
               offsetstart = graincentreradiusupdate-((graincentreradiusupdate**2)-((armwidthupdate/2)**2))**(1/2)
               startheight = armheightupdate+graincentreradiusupdate-offsetstart
               xstart = math.cos(thea)*(startheight)
               ystart = math.sin(thea)*(startheight)
               x.append(xstart), y.append(ystart)
              
               thea = ((armtrack - 1) * angle_between_arms) + math.pi/2
               startlength  = ((armwidthupdate)/2)-burn_fillet
               yrun = math.sin(thea)*(startlength)
               xrise = math.cos(thea)*(startlength)
               xendpoint = xstart+xrise
               yendpoint = ystart+yrun
               startarmx = numpy.linspace(xstart,xendpoint,100)
               startarmy = numpy.linspace(ystart,yendpoint,100) #breaks the arm into 100 points
               for xpoint in startarmx:
                   x.append(xpoint) 
               for ypoint in startarmy:
                   y.append(ypoint)
               xcurrent = xpoint
               ycurrent = ypoint
               ypoint2 =ypoint - armwidthupdate+burn_fillet*2
               
               #burn radius calculations
               if run == 0:
                   curve = numpy.linspace(0, burn_fillet ,500)
                   for curvey in curve:
                       curvex = ((burn_fillet**2)-(curvey**2))**(1/2)
                       ystore = ypoint+curvey
                       xstore = startheight-(burn_fillet-curvex)
                       r = ((xstore**2)+(ystore**2))**(1/2)
                       theac = math.acos(xstore/r)
                       radiuscurve.append(r)
                       theacurve.append(theac)
                       
                   curve2 = numpy.linspace(burn_fillet,0 ,500)
                   for curvey in curve2:
                       curvex = ((burn_fillet**2)-(curvey**2))**(1/2)
                       ystore = ypoint2-curvey
                       xstore = startheight-(burn_fillet-curvex)
                       r = ((xstore**2)+(ystore**2))**(1/2)
                       theac = math.asin(ystore/r)
                       radiuscurve2.append(r)
                       theacurve2.append(theac)
                       run = 1
                       
               curvey = 0
               #burn fillet when needed
               if burn_fillet > 0:
                   for i in range (0, 500):
                       angle_update = (armtrack - 1) * angle_between_arms
                       thea2 = theacurve[i] + angle_update
                       xcurrent = math.cos(thea2)*radiuscurve[i]
                       ycurrent = math.sin(thea2)*radiuscurve[i]
                       x.append(xcurrent)
                       y.append(ycurrent)
               #drop for arm        
               thea1 = (angle/2) 
               thea = (armtrack - 1) * angle_between_arms
               drop = armheightupdate
               yrun = math.sin(thea)*(drop)
               xrise = math.cos(thea)*(drop)
               xendpoint = xcurrent-xrise
               yendpoint = ycurrent-yrun
               startarmx = numpy.linspace(xcurrent,xendpoint,100)
               startarmy = numpy.linspace(ycurrent,yendpoint,100)
               for xpoint in startarmx:
                   x.append(xpoint) 
               for ypoint in startarmy:
                   y.append(ypoint)
               xcurrent = xpoint
               ycurrent = ypoint
               
               #curve of the bore
               thea1 =(angle/2)+(angle_between_arms*(armtrack-1))
               thea2 = thea1 + anglebore  #calculates the angle of the curve
               curvethea = numpy.linspace(thea1,thea2,500)
               for thea in curvethea:
                   xcurrent = math.cos(thea)*graincentreradiusupdate
                   ycurrent = math.sin(thea)*graincentreradiusupdate
                   x.append(xcurrent), y.append(ycurrent)
               #angled second arm
               
               thea = ((armtrack) * angle_between_arms) 
               rise = (armheightupdate-burn_fillet)
               yrun = math.sin(thea)*(rise)
               xrise = math.cos(thea)*(rise)
               xendpoint = xcurrent+xrise
               yendpoint = ycurrent+yrun
               startarmx = numpy.linspace(xcurrent,xendpoint,100)
               startarmy = numpy.linspace(ycurrent,yendpoint,100)
               for xpoint in startarmx:
                   x.append(xpoint) 
           
               for ypoint in startarmy:
                   y.append(ypoint)
               xcurrent = xpoint
               ycurrent = ypoint  
               if burn_fillet > 0:
                   for i in range (0, 500):
                       angle_update = (armtrack) * angle_between_arms
                       thea2 = theacurve2[i] + angle_update
                       xcurrent = math.cos(thea2)*radiuscurve2[i]
                       ycurrent = math.sin(thea2)*radiuscurve2[i]
                       x.append(xcurrent)
                       y.append(ycurrent)
               
               thea1 = (angle/2)
               thea = ((armtrack) * angle_between_arms)+ math.pi/2
               startlength  = ((armwidthupdate)/2)-burn_fillet
               yrun = math.sin(thea)*(startlength)
               xrise = math.cos(thea)*(startlength)
               xendpoint = xcurrent+xrise
               yendpoint = ycurrent+yrun
               startarmx = numpy.linspace(xcurrent,xendpoint,100)
               startarmy = numpy.linspace(ycurrent,yendpoint,100) #breaks the arm into 100 points
               for xpoint in startarmx:
                   x.append(xpoint) 
               for ypoint in startarmy:
                   y.append(ypoint)
                   # angle_current_arm = 
                  # xpoint = math.sin(angle_between_arms)*armheightupdate
                  # ypoint = math.cos(angle_between_arms)*armheightupdate
                  # gradient = ypoint/xpoint
                  # curvethea = numpy.linspace(thea1,thea2,500)
               #armtrack = armtrack+1
   return x,y,x1,y1
''' 
#could not get to work
y0 = numpy.zeros(5) # = [0,0]
y0[0] = numberofarms
y0[1] = grain_length
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