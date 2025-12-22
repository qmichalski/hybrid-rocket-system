# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 22:11:32 2025

@author: Guillaume
"""
import math
import numpy
import matplotlib.pyplot as graph

def graingeometry_finocyl_plotter_inital_shape (grainlength,armheight,armwidth,numberofarms,graincentreradius,chamber_outer_radius,):
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