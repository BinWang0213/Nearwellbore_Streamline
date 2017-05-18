#########################################################################                                                                      	
#		(C) 2017 Department of Petroleum Engineering,					#						
#		Univeristy of Louisiana at Lafayette, Lafayette, US.			#
#																		#
# This code is released under the terms of the BSD license, and thus	#
# free for commercial and research use. Feel free to use the code into	#
# your own project with a PROPER REFERENCE. 							#
#                                                                      	#
# A near wellbore streamline tracking code								#            	
# Author: Bin Wang														#	
# Email: binwang.0213@gmail.com											#	
# Reference: Wang, B., Feng, Y., Du, J., et al. (2017) An Embedded 		#
# Grid-Free Approach for Near Wellbore Streamline Simulation.			#
# doi:10.2118/SPE-182614-MS												#
#########################################################################

import numpy as np
import matplotlib.pyplot as plt


############################### Grid Subdivision Method #################################
def PointOnUnitSquare(NSL,Endpoint=False):
    '''Genetrating points around a unit square
    '''
    Pts=np.zeros((NSL,2))
    NSL_edge=int(NSL/4)
    if Endpoint==True:
        dx=np.linspace(0.0,1.0,NSL_edge,endpoint=True)
    if Endpoint==False:
        dx=np.linspace(0.0,1.0,NSL_edge,endpoint=True)

    dx_rev=dx[::-1]
    for i in range(4):
        for j in range(NSL_edge):
            if (i==0):
                    Pts[j+i*NSL_edge,0],Pts[j+i*NSL_edge,1]=dx[j],0
            if (i==1):
                    Pts[j+i*NSL_edge,0],Pts[j+i*NSL_edge,1]=1,dx[j]
            if (i==2):
                    Pts[j+i*NSL_edge,0],Pts[j+i*NSL_edge,1]=dx_rev[j],1
            if (i==3):
                    Pts[j+i*NSL_edge,0],Pts[j+i*NSL_edge,1]=0,dx_rev[j]
    return Pts

def PointOnUnitEdge(NSL,Endpoint=False):
    '''Genetrating points around a unit Edge
    '''
    Pts=np.zeros((NSL,2))
    if Endpoint==True:
        dx=np.linspace(0.0,1.0,NSL,endpoint=True)
    if Endpoint==False:
        dx=np.linspace(0.03,0.97,NSL,endpoint=True)

    for i in range(NSL):
            Pts[i,0],Pts[i,1]=0,dx[i]
    return Pts

def CalcDist(Pts0=(0,0),Pts1=(1,1)):
    '''Calculating distance of two points
    '''
    return np.sqrt((Pts1[0]-Pts0[0])**2+(Pts1[1]-Pts0[1])**2)


def RotateSL(SL,Single=0,origin=(0.0,0.0),angle=np.pi/2):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    http://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
    """
    #angle=self.RotateAngle
        
    NSL=len(SL)
    SL_new=SL.copy()
    ox,oy=origin
        
    if (Single==0):
        for i in range(NSL):
            for j in range(10): #1 SL have 10 nodes
                px,py=SL[i][j]
                qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
                qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
                SL_new[i][j]=qx,qy
    elif (Single==1):
            for i in range(10): #1 SL have 10 nodes
                px,py=SL[i]
                qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
                qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
                SL_new[i]=qx,qy
    return SL_new
    
def TranslateSL(SL,Single=0,new_origin=(1.0,1.0),origin=(0.0,0.0)):
    """
    Translate the SL by a given origin
    Single=0  Treatment multiple streamline
    Single=1  Treatment single streamline
    """
    #new_origin=self.NewOrigin
        
    NSL=len(SL)
    SL_new=SL.copy()
        
    ox_new,oy_new=new_origin
    ox,oy=origin
        
    if (Single==0):
        for i in range(NSL):
            for j in range(10): #1 SL have 10 nodes
                px,py=SL[i][j]
                qx = px +  (ox_new - ox)
                qy = py +  (oy_new - oy)
                SL_new[i][j]=qx,qy
        
    elif(Single==1):
        for i in range(10): #1 SL have 10 nodes
                px,py=SL[i]
                qx = px +  (ox_new - ox)
                qy = py +  (oy_new - oy)
                SL_new[i]=qx,qy
    return SL_new



############################### Embedded Method #################################
def cosspace(st,ed,N,endpoint=True):
    """
    Auto line segment refinement at end point
    e.g. --- - -  -  -  - - ---
    """
    #N=N+1
    AngleInc=np.pi/(N-1)
    CurAngle = AngleInc
    space=np.linspace(0,1,N,endpoint=endpoint)
    space[0]=st
    for i in range(N-1):
        space[i+1] = 0.5*np.abs(ed-st)*(1 - np.cos(CurAngle));
        CurAngle += AngleInc
    if ed<st:
        space[0]=ed
        space=space[::-1]

    return space

def centroid2D(Pts):
    Pts = np.asarray(Pts)
    Npts=len(Pts)
    return np.sum(Pts[:,0])/Npts,np.sum(Pts[:,1])/Npts

def EndPointOnLine(Pts_a=(0,0),Pts_b=(0,1),Nseg=4,refinement="linspace",Endpoint=True):
    '''Genetrating endpoints along a line segment
       algorithm: point with a given distance along a line: x=x_Start+unit_vector_x*distance  y=y_Start+unit_vector_y*distance
    Arguments
        ---------
        Pts_a    -- The start-point.
        Pts_b    -- The end-point
        Npts     -- Number of endpoints 
        Nseg     -- Number of segments
        unit_vx  -- Unit vector for x coordinates
        unit_vy  -- Unit vector for y coordinates
        interval -- Segment interval
                    uniform    - - - - - - - - - -  (linspace)
                    refinement -- - -  -  -  - - -- (cosspace)
    '''
    Npts=Nseg+1
    Pts=np.zeros((Npts,2))
    length=CalcDist(Pts_a,Pts_b)
    unit_vx=(Pts_b[0]-Pts_a[0])/length
    unit_vy=(Pts_b[1]-Pts_a[1])/length
    
    if (refinement=="linspace"):
        interval=np.linspace(0.0,length,Npts,endpoint=Endpoint)
        rinterval=np.linspace(length,0.0,Npts,endpoint=Endpoint)
    elif (refinement=="cosspace"):
        interval=cosspace(0.0,length,Npts,endpoint=Endpoint)
        rinterval=cosspace(length,0.0,Npts,endpoint=Endpoint)
    
    for i in range(Npts):
        Pts[i,0]=Pts_a[0]+interval[i]*unit_vx
        Pts[i,1]=Pts_a[1]+interval[i]*unit_vy
    
    return Pts

def EndPointOnCircle(Origin=(0,0),R=1,Nseg=4):
    '''Genetrating endpoints along a circle
    Arguments
        ---------
        Origin  -- The start-point.
        R       -- The end-point
        Npts    -- Number of endpoints 
        Nseg    -- Number of segments
    '''
    Npts=Nseg+1
    Pts=np.zeros((Npts,2))
    
    interval=np.linspace(0, 2*np.pi, Npts)

    for i in range(Npts):
        Pts[i,0]=Origin[0]+np.cos(interval[i])*R
        Pts[i,1]=Origin[1]+np.sin(interval[i])*R
    
    return Pts

def LineSegIntersect(Line1=([0,0],[1,1]),Line2=([0,1],[1,0])):
    #Algorithm from http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    #Test whether 2 line segment are intersected
    xa,ya,xb,yb=Line1[0][0],Line1[0][1],Line1[1][0],Line1[1][1]
    xc,yc,xd,yd=Line2[0][0],Line2[0][1],Line2[1][0],Line2[1][1]
    ccw_ACD=(yd-ya)*(xc-xa) > (yc-ya)*(xd-xa)
    ccw_BCD=(yd-yb)*(xc-xb) > (yc-yb)*(xd-xb)
    ccw_ABC=(yc-ya)*(xb-xa) > (yb-ya)*(xc-xa)
    ccw_ABD=(yd-ya)*(xb-xa) > (yb-ya)*(xd-xa)
    return ccw_ACD != ccw_BCD and ccw_ABC != ccw_ABD

############################### Fill-Grid Method #################################
def PolygonArea(Pts):
    #Calculate 2D polygon area using Shoelace formula
    #http://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates/30408825
    x,y=np.asarray(Pts)[:,0],np.asarray(Pts)[:,1]
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

