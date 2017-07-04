#########################################################################                                                                       
#       (C) 2017 Department of Petroleum Engineering,                   #                       
#       Univeristy of Louisiana at Lafayette, Lafayette, US.            #
#                                                                       #
# This code is released under the terms of the BSD license, and thus    #
# free for commercial and research use. Feel free to use the code into  #
# your own project with a PROPER REFERENCE.                             #
#                                                                       #
# A near wellbore streamline tracking code                              #               
# Author: Bin Wang                                                      #   
# Email: binwang.0213@gmail.com                                         #   
# Reference: Wang, B., Feng, Y., Du, J., et al. (2017) An Embedded      #
# Grid-Free Approach for Near Wellbore Streamline Simulation.           #
# doi:10.2118/SPE-182614-MS                                             #
#########################################################################

import numpy as np
import matplotlib.pyplot as plt
import Lib.smallestenclosingcircle as enclose #smallestenclosingcircle library from Nayuki https://www.nayuki.io/page/smallest-enclosing-circle
from .geometry import *

                            # Embedded Method #

######################
#
#  Basic Element Class
#
######################
class Panel:
    """Contains information related to a panel(Virtual Boundary elements)."""
    def __init__(self, xa, ya, xb, yb,Q_bd,P_bd,marker):
        """Creates a panel from A-B. It also can be used for real boundary with ghost pairs
        
        Arguments
        ---------
        xa, ya     -- Cartesian coordinates of the first end-point A.
        xb, yb     -- Cartesian coordinates of the second end-point B.
        xc, yc     -- Cartesian coordinates of the center point.
        length     -- length of this BE.
        sinalpha   -- sin(aj) in Fig. 1
        cosalpha   -- cos(aj) in Fig. 1
        Q_bd       -- Neumann BC the noraml pressure gradient
        P_bd       -- Dirchlet BC the pressure
        marker     -- boundary marker. e.g bd1,bd2,bd3...
        rab        -- half distance between two ghost pairs
        """
        self.xa, self.ya = xa, ya
        self.xb, self.yb = xb, yb
        
        self.xc, self.yc = (xa+xb)/2, (ya+yb)/2       # control-point (center-point)
        self.length = np.sqrt((xb-xa)**2+(yb-ya)**2)     # length of the panel
        
        #unit vector (x2-x1)/length  (y2-y1)/length
        #normal unit vector (y1-y2)/length   (x2-x1)/length
        #point with a given distance along a line: x?=x_Start+unit_vector_x*distance  y?=y_Start+unit_vector_y*distance
        self.rab=0.00003
        self.x_a=self.xc-self.rab/2*(ya-yb)/self.length
        self.y_a=self.yc-self.rab/2*(xb-xa)/self.length
        self.x_b=self.xc+self.rab/2*(ya-yb)/self.length
        self.y_b=self.yc+self.rab/2*(xb-xa)/self.length
        
        # orientation of the panel (angle between x-axis and panel)
        self.sinalpha=(yb-ya)/self.length
        self.cosalpha=(xb-xa)/self.length
        
        self.Qbd=Q_bd                              # Newmann boundary condition
        self.Pbd=P_bd                              # Dirchlet boundary condition
        self.Q = 0.                                # source strength
        self.P = 0.                                # source strength
        self.marker=marker                         # boundary marker


class Well:
    """Contains information related to a Well."""
    def __init__(self, xw, yw,rw,Q,P):
        """Creates a panel.
        
        Arguments
        ---------
        xw, yw -- Cartesian coordinates of well source.
        Q      -- Flow rate of well source/sink
        P      -- Bottom hole pressure of well source/sink
        rw     -- radius of well source.
        """
        self.xw, self.yw = xw, yw
        
        self.Q = Q                               # Newmann boundary condition
        self.P = P                               # Dirichelt boundary condition
        self.rw = rw                             # raidus

######################
#
#  Mesh Tools
#
######################

def Add_Line(Pts_a=(0,0),Pts_b=(0,1),Nbd=5,Qbd=10.0,Pbd=10.0,bd_marker=-1,panels=[]):
    """Creates a BE along a line boundary. Anticlock wise, it decides the outward normal direction
        
        Arguments
        ---------
        xa, ya -- Cartesian coordinates of the first start-point.
        xb, yb -- Cartesian coordinates of the second end-point
        Nbd    -- Number of elements in this line
        Qbd    -- Neumann B.C
        Pbd    -- Dirchelet B.C
        Panels -- BE container, return value
    """
    
    Pts=EndPointOnLine(Pts_a,Pts_b,Nseg=Nbd,refinement="cosspace")
    for i in range(Nbd):
        panels.append(Panel(Pts[i][0],Pts[i][1], Pts[i+1][0], Pts[i+1][1],Qbd,Pbd,bd_marker))#Neumann, Dirchlet

def Add_Circle(Pts_c=(0,0),R=1,Nbd=5,Qbd=10.0,Pbd=10.0,bd_marker=-1,panels=[]):
        """Creates a BE along a circle boundary. Clockwise for well discretization, it decides the outward normal direction
        
        Arguments
        ---------
        xc, yc -- Cartesian coordinates of the circle center.
        R      -- The radius of the circle
        Nbd    -- Number of elements in this line
        Qbd    -- Neumann B.C
        Pbd    -- Dirchelet B.C
        Panels -- BE container, return value
        """
    # Generalize the image well position/Define the image panel around the computation domain
        Pts=EndPointOnCircle(Pts_c,R=R,Nseg=Nbd)
        
        for i in range(Nbd):
            panels.append(Panel(Pts[i][0],Pts[i][1], Pts[i+1][0], Pts[i+1][1],Qbd,Pbd,bd_marker))

def Add_Well(loc=(0,0),rw=1,Qwell=0,Pwell=100,wells=[]):
        """Creates a BE along a circle boundary. Clockwise for well discretization, it decides the outward normal direction
        
        Arguments
        ---------
        loc   -- location of the well.
        rw    -- The radius of wellbore
        Qwell -- Wellbore flow rate
        Pwell -- Bottomhole pressure
        wells -- Well container, return value
        """
        wells.append(Well(loc[0],loc[1],rw,Qwell,Pwell))


######################
#
#  Core Code
#
######################
class WellGrid:
    """Contains information and functions related to a wellgrid."""
    def __init__(self, Pts_e=[],Pts_w=[],Qw=[],Qe=[],Nbd=10,rw=0.25,h=26.25,phi=0.2,miu=1,kxy=(100,100)):
        """Creates a WellGrid with abrituary shape and multiple wells using VBEM
        
        Arguments
        ---------
        Ne     -- Number of wellblock edges
        Nw     -- Number of wells
        Nbd    -- Number of virtual boundary elements per edge (ghost point pair).
        Nbe    -- Total number of boundary elements(unknows)
        Qw     -- Wellbore flow rate, stb/day
        Qe     -- Outflow edge flow rate, stb/day
        h      -- Wellblock thickness
        phi    -- Wellblock porosity
        miu    -- Fluid viscosity
        kxy    -- Principle permeability of the wellblock. e.g. [kxx,kyy]
        
        Pts_e  -- Vertices of boundary edge e.g. [(0,0),(0,1),(1,1),(1,0)]
        Pts_w  -- Location of well e.g. [(0.5,0.5),(0.75,0.75)]
        domain_min  -- minimum (x,y) for a polygon domain
        domain_max  -- maximum (x,y) for a polygon domain
        Pts Location Scheme example   
            1-----2 
            |     |             
            | 1.. |             
            |     |             
            0-----3             
        Pts        OutFlow Edge
        0          0(0,1)
        1          1(1,2)
        2          2(2,3)
        3          3(3,0)
        Well
        1
        2
        ..
        Default Sequence in wellgrid: bottom(0)-left(1)-top(2)-right(3), clockwise
        
        Derived&Output Arguments
        ---------
        BEs    -- List of circle virtual boundary elements [virtual boundary]
        Ghos   -- List of boundary ghost node pairs [real boundary]
        Wells  -- List of wells 

        SL     -- Array of streamlines
        TOF    -- Array of time-of-flight
        """
        self.Pts_e=Pts_e
        self.Pts_w=Pts_w
        self.Ne=len(Qe)
        self.Nw=len(Qw)
        self.Nbd=Nbd
        self.Nbe=self.Nbd*self.Ne
        
        self.Qw=Qw
        self.Qe=Qe
        self.rw=rw
        self.h=h
        self.phi=phi
        self.miu=miu
        self.kx=kxy[0]
        self.ky=kxy[1]
        
        #Boundary elements
        self.BEs=[]
        self.Ghos=[]
        self.Wells=[]
        
        #Streamlines
        self.SL=[]
        self.TOF=[]
        
        #Additional Geometory variable
        self.domain_min=(min(np.asarray(self.Pts_e)[:,0]),min(np.asarray(self.Pts_e)[:,1]))
        self.domain_max=(max(np.asarray(self.Pts_e)[:,0]),max(np.asarray(self.Pts_e)[:,1]))

    def Meshing(self):
        """Genetrating meshing for VBEM (virtual boundary elements, ghost node pairs and wells)
        Fig. 2 in in SPE-182614-MS
        
        Arguments
        ---------
        Origin     -- Centroid of a polygon domain
        R          -- Minimum radius of circle which enclose the polygon domain 
        error      -- empty space for plotting the Mesh
        """

        #Circle virtual boundary elements
        Origin = centroid2D(self.Pts_e)
        R = 1.5*enclose.make_circle(self.Pts_e)[1]
        #R=50
        #print(R)
        Add_Circle(Origin,R,self.Nbe,panels=self.BEs)
        
        #Ghost node pairs and Boundary conditions
        self.Pts_e.append(self.Pts_e[0])
        for i in range(self.Ne):
            Add_Line(self.Pts_e[i],self.Pts_e[i+1],self.Nbd,Qbd=self.Qe[i],panels=self.Ghos,bd_marker=i)
        
        #Wells
        for i in range(self.Nw):
            Add_Well(self.Pts_w[i],self.rw,Qwell=self.Qw[i],wells=self.Wells)
        
        #Mesh Plot
        error=R*0.2 #empty space around the circle VBEM

        plt.figure(figsize=(3, 3))
        plt.axes().set(xlim=[Origin[0]-R-error, Origin[0]+R+error],
                       ylim=[Origin[1]-R-error, Origin[1]+R+error],aspect='equal')

            #Domain boundary
        plt.plot(*np.asarray(self.Pts_e).T,lw=1,color='black')
        plt.scatter(*np.asarray(self.Pts_w).T,s=20,color='red')
        
            #Virtual Boundary elements
        plt.plot(np.append([BE.xa for BE in self.BEs], self.BEs[0].xa), 
             np.append([BE.ya for BE in self.BEs], self.BEs[0].ya), 
             'bo--',markersize=5);
        
            #Ghost node pair
        plt.scatter([Gho.x_a for Gho in self.Ghos], [Gho.y_a for Gho in self.Ghos], color='r',facecolors='none', s=5)
        plt.scatter([Gho.x_b for Gho in self.Ghos], [Gho.y_b for Gho in self.Ghos], color='r',facecolors='none', s=5)
            
        plt.axes().set_title('VBEM Mesh')
        plt.show()
    

    def GH_analytical(self,Pts=(0,0),BE=[]):
        '''Calclate BE influence coefficient for the pressure and normal flux
           Eq. (5) (6),(7),(8) in SPE-182614-MS
        '''
        #Transfer global coordinate point(x,y) to local coordinate
        x,y=Pts[0]-BE.xa,Pts[1]-BE.ya #Eq. A-1 in SPE-182614-MS
        L=BE.length
        kr=self.kx/self.ky
        unit_v=0.4468 #unit converstion factor

        a=BE.cosalpha**2+kr*BE.sinalpha**2
        b=x*BE.cosalpha+kr*BE.sinalpha*y
        c=y*BE.cosalpha-x*BE.sinalpha
        #dp=-70.6*self.miu/self.h/np.sqrt(self.kx*self.ky)
        dv=-unit_v/self.h/self.phi*np.sqrt(kr)  
        
        #print('xy',x,y)
        #print('abc',a,b,c)
        #print('kr,sin,cos',kr,BE.cosalpha,BE.sinalpha)
        Gij = -1/a*(
                             (
                              b*np.log(x**2-2*b*L+a*L**2+kr*y**2)
                             -L*a*np.log((x-L*BE.cosalpha)**2+kr*(y-L*BE.sinalpha)**2)
                             +2*np.sqrt(kr)*c*np.arctan((b-a*L)/np.sqrt(kr)/c)
                             )
                             -
                             (
                               b*np.log(x**2+kr*y**2)
                               +2*np.sqrt(kr)*c*np.arctan((b)/np.sqrt(kr)/c)
                             )         
                )
        Hij_x = dv/a*(
                             ( 
                            BE.cosalpha*np.log(x**2-2*b*L+a*L**2+kr*y**2)
                            + 2*np.sqrt(kr)*BE.sinalpha*np.arctan((a*L-b)/np.sqrt(kr)/c) 
                             )
                            -
                             (
                            BE.cosalpha*np.log(x**2+kr*y**2)+2*np.sqrt(kr)*BE.sinalpha*np.arctan((-b)/np.sqrt(kr)/c)
                             )    
                     )  
        Hij_y = dv/a*(
                           ( 
                            BE.sinalpha*np.log(x**2-2*b*L+a*L**2+kr*y**2)
                            + 2*np.sqrt(1/kr)*BE.cosalpha*np.arctan((b-a*L)/np.sqrt(kr)/c) 
                             )
                            -
                             (
                            BE.sinalpha*np.log(x**2+kr*y**2)+2*np.sqrt(1/kr)*BE.cosalpha*np.arctan((b)/np.sqrt(kr)/c)
                             )    
                     ) 
        return Gij,Hij_x,Hij_y
    
    def GHw_analytical(self,Pts=(0,0),well=[]):
        '''Calclate Well influence coefficient for the pressure and normal flux
           Eq. (9),(10),(11) in SPE-182614-MS
        '''
        unit_v=0.8936 #0.8936 is unit converstion factor

        #dp=-70.6*self.miu/self.h/np.sqrt(self.kx*self.ky)
        dv=unit_v/self.h/self.phi*np.sqrt(self.kx/self.ky)
        
        Gij=np.log((Pts[0]-well.xw)**2+(self.kx/self.ky)*(Pts[1]-well.yw)**2)
        Hij_x=dv*(Pts[0]-well.xw)/((Pts[0]-well.xw)**2+(self.kx/self.ky)*(Pts[1]-well.yw)**2)
        Hij_y=dv*(Pts[1]-well.yw)/((Pts[0]-well.xw)**2+(self.kx/self.ky)*(Pts[1]-well.yw)**2)
        return Gij,Hij_x,Hij_y
    
    def FlowSol(self):
        '''Solve the flow field (build matrix) using VBEM
           Eq. (13) in SPE-182614-MS
        '''
        
        #1. Build matrix A
        MatA=np.zeros((self.Nbe,self.Nbe))
        for i, Gho in enumerate(self.Ghos): #target ghost pairs nodes
            for j, BE in enumerate(self.BEs): #BE source
                Pts_a=Gho.x_a,Gho.y_a
                Pts_b=Gho.x_b,Gho.y_b
                MatA[i,j]=self.GH_analytical(Pts_a,BE)[0]-self.GH_analytical(Pts_b,BE)[0]
        
        #2. Build matrix RHS
        RHS_well=np.zeros((self.Nbe))
        RHS_Neumann=np.zeros((self.Nbe))
        
        for i, Gho in enumerate(self.Ghos): #target ghost pairs nodes
            tempRHS=0.0
            for j, Well in enumerate(self.Wells): #BE source
                Pts_a=Gho.x_a,Gho.y_a
                Pts_b=Gho.x_b,Gho.y_b
                tempRHS=tempRHS+(self.GHw_analytical(Pts_a,Well)[0]-self.GHw_analytical(Pts_b,Well)[0])*Well.Q
            RHS_well[i]=tempRHS
        for i in range(self.Ne):#Boundary conditions
            for j, Gho in enumerate(self.Ghos):#Corrsponding elements
                if(Gho.marker==i):
                    #print(i,Gho.Qbd)
                    Lbd=CalcDist(self.Pts_e[i],self.Pts_e[i+1])
                    kbd=self.kx*abs(Gho.sinalpha)+self.ky*abs(Gho.cosalpha)
                    if (self.kx==self.ky): kbd=self.kx
                    RHS_Neumann[j]=4*np.pi*Gho.Qbd*Gho.rab*np.sqrt(self.kx*self.ky)/Lbd/kbd
        RHS=-RHS_well-RHS_Neumann
        
        #4. Solve the matrix
        Q_BEs=np.linalg.lstsq(MatA,RHS)[0] #Solve Equation
        for i, BE in enumerate(self.BEs):
            BE.Q=Q_BEs[i]
        
        #np.savetxt('MatA.csv',MatA)
        #np.savetxt('Matb.csv',RHS)
        #print(Q_BEs)
    

    def FieldSol(self,Pts=(0,0)):
        '''Calculate the pressure and velocity at any point (x,y)
           Eq. (2) (3) and (4) in SPE-182614-MS
        '''
        unit_p=70.6 #unit converstion factor

        dp=-unit_p*self.miu/self.h/np.sqrt(self.kx*self.ky)
        p=u=v=0.0
        for i,BE in enumerate(self.BEs):
                puv=self.GH_analytical(Pts,BE)
                p=p+BE.Q*dp*puv[0]
                u=u+BE.Q*puv[1]
                v=v+BE.Q*puv[2]
        for i,Well in enumerate(self.Wells):
                puv_w=self.GHw_analytical(Pts,Well)
                p=p+Well.Q*dp*puv_w[0]
                u=u+Well.Q*puv_w[1]
                v=v+Well.Q*puv_w[2]
        return p,u,v

    def FieldPlot(self,vmax=100.01):
        '''Plot pressure&velocity field and Preview the streamline
        '''

        #Calculate pressure and velocity field
        from matplotlib import path
        Polygon = path.Path(self.Pts_e)
        
        N = 30                  # number of points in the x and y directions
        error=1e-6
        xmin,ymin=self.domain_min[0]+error,self.domain_min[1]+error
        xmax,ymax=self.domain_max[0]-error,self.domain_max[1]-error
        X, Y = np.meshgrid(np.linspace(xmin, xmax, N), np.linspace(ymin, ymax, N))  # generates a mesh grid
        #Calculate the velocity and pressure field
        p = np.empty((N, N), dtype=float)
        u = np.empty((N, N), dtype=float)
        v = np.empty((N, N), dtype=float)
        
        for i in range(N):
            for j in range(N):
                Pts=(X[i,j], Y[i,j])
                flag=Polygon.contains_points([Pts])
                #print(Pts,flag)
                #flag=True
                if (flag==True):
                    puv=self.FieldSol(Pts)
                    p[i,j],u[i,j],v[i,j]=puv[0],puv[1],puv[2]
                else:#point is not within the domain
                    p[i,j]=u[i,j]=v[i,j]= "nan"
        

        
        fig, axes = plt.subplots(ncols=3,figsize=(10, 10))
        Vtotal= np.sqrt(u**2+v**2)
        if (vmax==100.01):
            Vtotal_max=np.nanmax(Vtotal.flatten())
        else:
            Vtotal_max=vmax
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        for i, ax in enumerate(axes.flat):
            ax.set(xlim=[xmin-error, xmax+error],ylim=[ymin-error, ymax+error],aspect='equal')
            ax.plot(*np.asarray(self.Pts_e).T,lw=1,color='black')
            ax.scatter(*np.asarray(self.Pts_w).T,s=20,color='red')
            if i==0:
                ax.set_title(r'Velocity Contour')
                level = np.linspace(0, Vtotal_max, 15, endpoint=True)
                im=ax.contour(X, Y, Vtotal,level,linewidths=1.2)
            if i==1:
                ax.set_title(r'Velocity Field')
                #im=ax.pcolormesh(X,Y,Vtotal,vmax=Vtotal_max)
                extent=(xmin,xmax,ymin,ymax)
                im=ax.imshow(Vtotal,vmin=0,vmax=Vtotal_max,extent=extent,origin='lower',interpolation='nearest')
            if i==2:
                import matplotlib.colors as colors
                import warnings
                warnings.filterwarnings('ignore') #hide the warnning when "nan" involves
                ax.set_title(r'Streamline Preview')
                level=colors.Normalize(vmin=0,vmax=Vtotal_max)
                strm=ax.streamplot(X, Y, u, v,color=Vtotal,norm=level)
                im=strm.lines
            #Add colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", "10%", pad=0.15)
            fig.colorbar(im,cax=cax)  # draw colorbar
        
        fig.tight_layout()
        plt.savefig('Field Plot.png',dpi=300)
        plt.show()
        return p,u,v
    
    def RungeKutta(self,Pts1=(0,0),TOF1=0.0,dt=0.0003,method="RK2",tol=0.05,debug=0):
        """Runge-kutta method for tracing streamline
           Method list:
           1. 2th order Runge-kutta method
           2. 4th order Runge-kutta method
           3. Adaptive method with given maximum distance for a specific dt
              maximum travel distance=tol*the minimum radius of enclose circle
        """
        
        
        if (method=="RK2"):
            #Step1
            puv1=self.FieldSol(Pts1)
            Vx1=[puv1[1],puv1[2]]
            Pts2=Pts1+np.multiply(Vx1,dt)
            #print(Pts1,Vx1,dt,Pts2)
            #Step2
            puv2=self.FieldSol(Pts2)
            Vx2=[puv2[1],puv2[2]]
            Pts_RK2=Pts1+np.multiply(np.add(Vx1,Vx2),0.5*dt)
            #print(Pts1,Vx2,dt,Pts_RK2)
            Pts=Pts_RK2
            
        
        if (method=="RK4"):#Eq. C1-C6
            #Step1
            puv1=self.FieldSol(Pts1)
            Vx1=[puv1[1],puv1[2]]
            kuv1=Vx1*dt
            Pts2=Pts1+0.5*kuv1
            
            #Step2
            puv2=self.FieldSol(Pts2)
            Vx2=[puv2[1],puv2[2]]
            kuv2=Vx2*dt
            Pts3=Pts1+0.5*kuv2
            
            #Step3
            puv3=self.FieldSol(Pts3)
            Vx3=[puv3[1],puv3[2]]
            kuv3=Vx3*dt
            Pts4=Pts1+0.5*kuv3
            
            #Step4
            puv4=self.FieldSol(Pts4)
            Vx4=[puv4[1],puv4[2]]
            kuv4=Vx4*dt
            
            Pts_RK4=Pts1+(kuv1+2*kuv2+2*kuv3+kuv4)/6
            Pts=Pts_RK4
        
        if (method=='Adaptive'):
            maxdist=enclose.make_circle(self.Pts_e)[1]*tol
            t=0.0
            count=0
            dt_sub=0.0
            #print('Start Point',Pts1)
            while (t<dt): #tracing streamline until it hit the boundary
                puv=self.FieldSol(Pts1)
                Vx1=[puv[1],puv[2]]
                
                V_Pts=np.sqrt(puv[1]**2+puv[2]**2)
                dt_sub=maxdist/V_Pts
                Pts_adaptive=Pts1+np.multiply(Vx1,dt_sub)
                #print('from',Pts1,'to',Pts_adaptive,'V',V_Pts,'dt',dt_sub,'t',t)
                t=t+dt_sub
                count=count+1
                Pts1=Pts_adaptive
                if (count>30):
                    #print('Diverging!')
                    break
            #print('Velocity',V_Pts,'MiniR',maxdist/tol,'MaximumDist',maxdist)
            if (debug): print('EndPoint',Pts_adaptive,'sub-dt',dt_sub,'sub-steps',count)
            Pts=Pts_adaptive
            dt=t
        
        
        TOF=np.add(TOF1,dt)
            
        return Pts,TOF
            
    def SLtrace(self,NSL=10,deltaT=0.1,method='adaptive',tol=0.05,debug=0):
        """Trace Streamlines in the wellgrid using runge-kutta method
        Arguments
        ---------
        NSL        -- Total Number of Streamline
        NSL_w      -- Streamline number for each well
        method     -- Numerical intergation method. RK2, RK4 and Adaptive are provided

        Output
        ---------
        TOF        -- TOF array of of the [NSL] streamlines  
        SL         -- nodes array of of the [NSL] streamlines  
        TOF_end    -- a boundary TOF list of the [NSL] streamlines  
        SL_end     -- a boundary node list of the [NSL] streamlines 
        """
        self.SL=[]
        self.TOF=[]
        TOF_end=[]
        SL_end=[]
        #Genetrating start points of streamline
        NSL_w=np.zeros((self.Nw))
        SL_count=0
        for i in range(self.Nw):
            NSL_w[i]=int(self.Qw[i]/sum(self.Qw)*NSL)
            for j in range(int(NSL_w[i])):
                Pts0=[[self.Wells[i].xw+self.Wells[i].rw*np.cos(j*2*np.pi/NSL_w[i]),
                       self.Wells[i].yw+self.Wells[i].rw*np.sin(j*2*np.pi/NSL_w[i])]]
                TOF0=[[0.0]]
                self.SL.append(Pts0)
                self.TOF.append(TOF0)
        
        NSL=sum(NSL_w) #update the streamline number
        
      
        #Tracing SL using RK algorithm
        for i in range(int(NSL)):
            j=0
            flag=False
            Pts0=(0,0)
            TOF0=0.0
            while (flag==False): #tracing streamline until it hit the boundary
                Pts0=self.SL[i][j]
                TOF0=self.TOF[i][j]
                
                Pts1,TOF1=self.RungeKutta(Pts0,TOF0,dt=deltaT,method=method,tol=tol,debug=debug)
                self.SL[i].append(Pts1)
                self.TOF[i].append(TOF1)
            
                j=j+1
                #print('j',j,self.SL)
                #print(self.TOF)
                #Check boundary hit
                for k in range(self.Ne):
                    SLseg=(Pts0,Pts1)
                    BD=(self.Pts_e[k],self.Pts_e[k+1])
                    flag=LineSegIntersect(SLseg,BD)
                    if (flag==True): break
                if (flag==True):
                    TOF_end.append(TOF0)
                    SL_end.append(Pts0)
                #if j==10: break
        
        #Plot Streamline
        plt.figure(figsize=(3, 3))
        plt.axes().set(xlim=[self.domain_min[0], self.domain_max[0]],ylim=[self.domain_min[1], self.domain_max[1]],aspect='equal')
        plt.title('Streamlines')
        #Domain boundary
        plt.plot(*np.asarray(self.Pts_e).T,lw=1,color='black')
        plt.scatter(*np.asarray(self.Pts_w).T,s=20,color='red')
        #Streamline
        for i in range(len(self.SL)):
            plt.plot(*np.asarray(self.SL[i]).T,lw=1,marker='o',markersize=0,color='blue')
        
        plt.show()
        return self.SL,self.TOF,SL_end,TOF_end