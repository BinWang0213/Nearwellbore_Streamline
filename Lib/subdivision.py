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
from .geometry import *

                    # Grid Subdivision Method #
# This code is adopted from Chapter 8 of Streamline Simulation Theory 
# and Practice (Datta-Gupta and King, 2007)



###############################
#
#  Basic Traingular Grid
#
###############################
class Subgrid:
    """Contains information and functions related to a transformed Triangle Grid.
    """
    def __init__(self, theta=0, r=1,h=10,phi=0.2,Qa=(10,10),Qb=(10,10)):
        """Creates a Subgrid which transformed a 2D triangle grid to a unit square
        
        Arguments
        ---------
        theta  -- The degree of the angle located at the well, rad
        r      -- The radius of the sub-trigular-grid, ft
        Qa     -- Flow rate in the alpha direction
        Qa0    -- Wellbore flow rate, stb/day
        Qa1    -- Outflow boundary flow rate, stb/day
        Qb     -- Flow rate in the beta direction 
        Qb0    -- Flow rate on internal edge where beta=0 ,stb/day
        Qb1    -- Flow rate on internal edge where beta=1 ,stb/day
        Jacobian -- Unit=ft3
        c1     -- Flux gradient in alpha direction
        c2     -- Flux gradient in beta direction
        
        Subscripts
        ---------
        1,a    -- alpha direction (x direction)
        2,b    -- beta direction (y direction)
        
        NewOrigin   --Auxilary parameters for Translation
        RotateAngle --Auxilary paramters for Rotation
        """
        self.theta=theta
        self.r=r
        self.phi=phi
        self.h=h
        
        self.Qa0=Qa[0]
        self.Qa1=Qa[1]
        self.Qb0=Qb[0]
        self.Qb1=Qb[1]
        
        self.c1=self.Qa1-self.Qa0
        self.c2=self.Qb1-self.Qb0
        
        self.NewOrigin=(25,25)
        self.RotateAngle=-np.pi*3/4
        
        self.Te=0
    
    def CalcT(self,c,inte_pair=(0,1)):
        """Calculate Pusedo time-of-flight in transformed space
           Eq. 8.25 in Page 251 (Datta-Gupta and King, 2007)
        """

        Q0=inte_pair[0] #Flux of particle-Lower limit of integator
        Q1=inte_pair[1] #Flux of target end face - Upper limit of integator
        
        if (Q1!=0 and Q0!=0): 
            T=1/c*np.log((Q1)/(Q0))
        else: #Case of no flow boundary 
            T=9e9

        return T
    
    def CalcPts(self,T,Pts=(0,0)):
        """Calculate particle location in transformed space
           Eq. 8.27 in Page 251 (Datta-Gupta and King, 2007)
        """

        Qa,Qb=self.LinearFlux(Pts[0],Pts[1])
                
        if (self.c1!=0):
            End_x=Pts[0]+Qa*(1/self.c1)*(np.exp(self.c1*T)-1)
        else: #Case 3(4) in SLTrace
            End_x=Pts[0]+Qa*T
        if (self.c2!=0):
            End_y=Pts[1]+Qb*(1/self.c2)*(np.exp(self.c2*T)-1)
        else: #Case 3(4) in SLTrace
            End_y=Pts[1]+Qb*T
            
        return End_x,End_y
    
    def Pts2Physic(self,Pts=(0,0)):
        """Transform the coordinates of a point from transfromed space to physical space
           Eq. 8.46 or 8.3 in Page 290,256 which vary with Jacobian
        """
        
        alpha,beta=Pts[0],Pts[1]
        x=alpha*self.r*(1-beta+beta*np.cos(self.theta))
        y=alpha*beta*self.r*np.sin(self.theta)
        return (x,y)
    
    def T2Physic(self,T=1,Pts=(0,0)):
        """Transform the pusedo time-of-flight(T) to physical time-of-flight(tau)
           Eq. 8.28 in Page 261 which vary with Jacobian. Analytical or numercial intergation required
        """
        alpha=Pts[0]
        Qa=self.LinearFlux(alpha)[0]
        
        if (self.c1!=0):
            tau_T=self.phi*(self.r**2)*self.h*np.sin(self.theta)*(alpha*T-Qa/self.c1*T+Qa/(self.c1**2)*np.exp(self.c1*T))
            tau_0=self.phi*(self.r**2)*self.h*np.sin(self.theta)*(Qa/(self.c1**2))
            tau=(tau_T-tau_0)*0.1781
        else:
            tau=self.phi*(self.r**2)*self.h*np.sin(self.theta)*(alpha*T+Qa/2*T*T)*0.1781
        
        return tau
        
    def LinearFlux(self,alpha=0,beta=0):
        """Calculate the linearly approximated flux
           Eq. 8.11 in Page 257
        """
        return self.Qa0+alpha*self.c1,self.Qb0+beta*self.c2
    
    def Trace1SL(self,Pts=(0,0.5),TOF_base=0.0,debug=0):
        """Trace a Streamline in the transformed space using Pollock's Algorithm
        
        Arguments
        ---------
        Pts      --  Initial particle point 
        T        --  The pesudo time-of-flight to reach four faces
        TOF_base --  Base Time-of-flight for connected streamlines, default=0 
        
        Subscripts
        ---------
        1,a    -- alpha direction (x direction)
        2,b    -- beta direction (y direction)

        T    0       1        2       3
           alpha0  alpha1   beta0   beta1
           well    outflow  bottom   top
        
        Normal Case
        ->1| |->2     ->1 != ->2
        
        Special Cases
        Case1 -Case(a) in Pollock's paper 
              -All inflow in a specific direction             -> | | <-
        Case2 -Case(b) in Pollock's paper
              -All outflow in a specific direction            <- | | -> 
        Case3 -Case(c) in Pollock's paper
              -Constant flow (velocity) in a specific direction      ->1| | ->2   ->1 = ->2
        Case4 -No flow in both faces in a specific direction    0| |0 
        """
        
        #---------------Flow case determination------------------
        CaseID_a=0 #Normal Case
        CaseID_b=0 #Normal Case
        if (self.Qa0>0.00001 and self.Qa1<-0.00001): CaseID_a=1
        if (self.Qb0>0.00001 and self.Qb1<-0.00001): CaseID_b=1
        if (self.Qa0<-0.00001 and self.Qa1>0.00001): CaseID_a=2
        if (self.Qb0<-0.00001 and self.Qb1>0.00001): CaseID_b=2
        if (abs(self.c1)<0.00001): CaseID_a=3
        if (abs(self.c2)<0.00001): CaseID_b=3
        if (abs(self.Qa0)<0.00001 and abs(self.Qa1)<0.00001): CaseID_a=4
        if (abs(self.Qb0)<0.00001 and abs(self.Qb1)<0.00001): CaseID_b=4
        
        T=np.ones(4)*123e9 
        Qa,Qb=self.LinearFlux(Pts[0],Pts[1])
        #---------------Streamline Tracing------------------
        #1. Pesudo Time of flight T in 4 direction
        if (CaseID_a==0):
            T[0]=self.CalcT(self.c1,(Qa,self.Qa0))
            T[1]=self.CalcT(self.c1,(Qa,self.Qa1))
        if (CaseID_a==1):
            T[0]=T[1]=1.1e9 #Particle never leave in this direction
        if (CaseID_a==2):
            if(Qa>0.0001): #if Qa>0 leave on the right face
                T[0]=2.1e9
                T[1]=self.CalcT(self.c1,(Qa,self.Qa1))
            if(Qa<-0.0001): #if Qa<0 leave on the left face
                T[0]=self.CalcT(self.c1,(Qa,self.Qa0))
                T[1]=2.2e9
            if (abs(Qa)<0.00001): #particle in stagenent region, never leave in this direction
                T[0]=T[1]=0
        if (CaseID_a==3):
            T[0]=(0-Pts[0])/Qa
            T[1]=(1-Pts[0])/Qa
        if (CaseID_a==4):
            T[0]=T[1]=4.1e9
        
        if (CaseID_b==0):
            T[2]=self.CalcT(self.c2,(Qb,self.Qb0))
            T[3]=self.CalcT(self.c2,(Qb,self.Qb1))
        if (CaseID_b==1):
            T[2]=T[3]=1.1e9 #Particle never leave in this direction
        if (CaseID_b==2):
            if(Qb>0.0001): #if Qb>0 leave on the upper face
                T[2]=2.1e9
                T[3]=self.CalcT(self.c2,(Qb,self.Qb1))
            if(Qb<-0.0001): #if Qb<0 leave on the lower face
                T[2]=self.CalcT(self.c2,(Qb,self.Qb0))
                T[3]=2.2e9
            if (abs(Qb)<0.00001): #particle in stagenent region, never leave in this direction
                T[2]=T[3]=0
        if (CaseID_b==3):
            T[2]=(0-Pts[1])/Qb
            T[3]=(1-Pts[1])/Qb
        if (CaseID_b==4):
            T[2]=T[3]=4.1e9
        
        #Special treatment of dt<=0
        if (Pts[0]==0 and Qa>0.00001): T[0]=9.1e9
        if (Pts[0]==1 and Qa<-0.00001): T[1]=9.1e9
        if (Pts[1]==0 and Qb>0.00001): T[2]=9.1e9
        if (Pts[1]==1 and Qb<-0.00001): T[3]=9.1e9
        
        #Special treatment of intial point located on the no-flow boundary
        if (Pts[0]==0 and Qa==0): T[0]=0.0
        if (Pts[0]==1 and Qa==0): T[1]=0.0
        if (Pts[1]==0 and Qb==0): T[2]=0.0
        if (Pts[1]==1 and Qb==0): T[3]=0.0

        
        for i in range(len(T)):
            if T[i]<-0.000001:
                T[i]=9e9
        
        if (debug): print("Flow Case (x,y)-(alpha,beta)",CaseID_a,CaseID_b)
        if (debug): print("Initial Point",Pts)
        
        
        Te=min(T)
        self.Te=Te

        if (debug): print(T)

        #2. Calculate End point 
        Pts_end=[0,0]
        Pts_end[0],Pts_end[1]=self.CalcPts(Te,Pts)
        
        if (debug): print("End Point",Pts_end,'End Time',self.Te)
        if (debug): print(T)
        
        #3.Calculate the streamline nodes and TOF in both coordinates
        SL_T=np.linspace(0,self.Te,10)
        SL_tau=np.zeros(len(SL_T))+TOF_base
        SL_phy=np.zeros((len(SL_T),2))
        SL_trans=np.zeros((len(SL_T),2))
                
        for i in range(len(SL_T)):
            SL_trans[i,0],SL_trans[i,1]=self.CalcPts(SL_T[i],Pts)
            SL_phy[i,0],SL_phy[i,1]=self.Pts2Physic(Pts=(SL_trans[i,0],SL_trans[i,1]))
            SL_tau[i]+=self.T2Physic(SL_T[i],Pts)
            
            #fix the tiny error, the number of e-17
            if (SL_trans[i,0]<0.0000001): SL_trans[i,0]=0.0
            if (SL_trans[i,1]<0.0000001): SL_trans[i,1]=0.0
            if (SL_phy[i,0]<0.0000001): SL_phy[i,0]=0.0
            if (SL_phy[i,1]<0.0000001): SL_phy[i,1]=0.0
        
        #SL_tau=SL_tau-SL_tau[0]
        
        return self.Te,Pts_end,SL_trans,SL_phy,SL_T,SL_tau
    
    def RotateSL(self,SL,Single=0,origin=(0.0,0.0),angle=np.pi/2):
        """
        Rotate streamlines counterclockwise by a given angle around a given origin (optional).
        
        The angle should be given in radians.
        http://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
        """
        angle=self.RotateAngle
        
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
    
    def TranslateSL(self,SL,Single=0,new_origin=(1.0,1.0),origin=(0.0,0.0)):
        """
        Translate the SL by a given origin
        Single=0  Treatment multiple streamline
        Single=1  Treatment single streamline
        """
        new_origin=self.NewOrigin
        
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
        
        
    def SLTrace(self,NSL=12,Pts=[]):
        """Construct and Plot streamline in both physical and transformed space
        
        Arguments
        ---------
        Pts - Initial particle start point (optional1)
        NSL - Number of Streamline (optional2)
        One has to made a choice of the two options

          Output
        ---------
        TOF        -- TOF array of of the [NSL] streamlines  
        SL         -- nodes array of of the [NSL] streamlines  
        TOF_end    -- a boundary TOF list of the [NSL] streamlines  
        SL_end     -- a boundary node list of the [NSL] streamlines 
        """
       
        #Grid edge
        Bound_vert=[(0,0),(1,0),(1,1),(0,1),(0,0)]
        Bound_vert_phy=[]
        for i in range(len(Bound_vert)):
            Bound_vert_phy.append(self.Pts2Physic(Bound_vert[i]))
        
        #Streamline
        if(len(Pts)==0): #if the initial Pts are not provided
            Pts=PointOnUnitSquare(NSL,Endpoint=False)
        else:
            NSL=len(Pts)
        
        SL=[]
        SL_phy=[]
        TOF_phy=[]
        
        for i in range(len(Pts)):
            temp=self.Trace1SL(Pts[i])
            SL.append(temp[2])
            SL_phy.append(temp[3])
            TOF_phy.append(temp[5])
        
        #SL_phy=self.RotateSL(SL_phy)
        #SL_phy=self.TranslateSL(SL_phy)
        
        fig, axs = plt.subplots(ncols=2)
        
        ax=axs[0]
        ax.plot(*np.asarray(Bound_vert).T,lw=3,color='red')
        for i in range(len(Pts)):
            ax.plot(*np.asarray(SL[i]).T,lw=1,marker='o',markersize=1,color='blue')
        ax.set_ylim(bottom=0)
        ax.set_aspect('equal')
        ax.set_title(r'Transformed Space ($\alpha,\beta$)')
        
        ax=axs[1]
        ax.plot(*np.asarray(Bound_vert_phy).T,lw=3,color='red')
        for i in range(len(Pts)):
            ax.plot(*np.asarray(SL_phy[i]).T,lw=1,marker='o',markersize=1,color='blue')
        ax.set_ylim(bottom=0)
        ax.set_aspect('equal')
        ax.set_title(r'Physical Space ($x,y$)')

        fig.tight_layout()
        plt.show()
        return SL_phy,TOF_phy


###############################
#
#  Wellblock with subdivision
#
###############################

class WellGrid:
    """Contains information related to a Triangle Grid."""
    def __init__(self, Rect0=(0,0),Rect1=(1,1),Qw=1000,Qe=(250,250,250,250),h=10,phi=0.2):
        """Creates a WellGrid with a center well (Gupta and King, 2010, Chapter 8)
        
        Arguments
        ---------
        Rect0  -- The coordinates of left bottom corner point of the rectangle
        Rect1  -- The coordinates of right top corner point of the rectangle
        Qw     -- Wellbore flow rate, stb/day
        Qe     -- Outflow edge flow rate, stb/day
        Qe_int -- Internal edge flow, stb/day
        h      -- Wellblock thickness
        phi    -- Wellblock porosity
        
        Pts Location Scheme   SubGrid ID Scheme     
            1-----2             .------. 
            |     |             |   2  |
            |  4  |             | 1 . 3| 
            |     |             |   0  |
            0-----3             .------.
        Pts        SubGrid     Neighbor    Internal Edge    OutFlow Edge
        0          0(3,4,0)     (1,3)        0(3,4)            0(3,0)
        1          1(0,4,1)     (2,0)        1(4.0)            1(0,1)
        2          2(1,4.2)     (3,1)        2(1,4)            2(1,2)
        3          3(2,4,3)     (0,2)        3(2,4)            3(2,3)
        4
        
        Default Sequence in both subgrid and wellgrid: bottom(0)-left(1)-top(2)-right(3), clockwise
        """
        #Cornor Points
        self.Pts=[(Rect0[0],Rect0[1]),(Rect0[0],Rect1[1]),(Rect1[0],Rect1[1]),(Rect1[0],Rect0[1])]
        
        #Center Well
        self.Pts.append(((Rect0[0]+Rect1[0])/2,(Rect0[1]+Rect1[1])/2) )
        
        self.Qw=Qw
        self.Qe=Qe
        self.Qe_int=np.zeros(4)
        self.h=h
        self.phi=phi
        self.theta=[np.pi/2,np.pi/2,np.pi/2,np.pi/2]
        
        #Subgrid
        self.SubGrids=[]
        self.NeighborID=[(1,3),(2,0),(3,1),(0,2)]
        
        #Streamline
        self.SL=[]
        self.TOF=[]
    
    def Subdivision(self,debug=0):
        '''Solving the flux balance equation and generating the subgrids (4 in this case)
           Eq. 8.54 in Page 292
        '''
        #1. Solving the flux balance equation
        MatQ=np.zeros((4,4))#Left hand Matrix
        np.fill_diagonal(MatQ,1)
        MatQ[0,1]=MatQ[1,2]=MatQ[2,3]=MatQ[3,0]=-1.0

        RHSQ=np.zeros(4) #Right hand side Qe-Qw
        for i in range(4):
            RHSQ[i]=self.Qe[i]-self.Qw*self.theta[i]/np.pi/2
        
        self.Qe_int=np.linalg.lstsq(MatQ,RHSQ)[0] #Solve Equation
        for i in range(4):
            if (abs(self.Qe_int[i])<0.00001):
                self.Qe_int[i]=0.0
        
        if (debug): print(MatQ,RHSQ,self.Qe_int)
        
        #2. Genetrating Subgrids by using calculated fluxes 
        theta0=self.theta[0]
        r0=CalcDist(self.Pts[3],self.Pts[4])
        self.Qe_int=np.append(self.Qe_int,self.Qe_int[0]) #for natural i+1 iteration
        
        for i in range(4):
            Qa=(self.Qw*theta0/2/np.pi,self.Qe[i])
            Qb=(-self.Qe_int[i+1],-self.Qe_int[i])
            if (debug): print(theta0,r0,Qa,Qb)
            self.SubGrids.append(Subgrid(theta=theta0,r=r0,h=self.h,phi=self.phi,Qa=Qa,Qb=Qb))
            if (debug): self.SubGrids[3].SLTrace(NSL=100)
            
        #3. Setup the grid oritation and origin
        self.SubGrids[0].NewOrigin=(25,25)
        self.SubGrids[0].RotateAngle=-np.pi*3/4
        self.SubGrids[1].NewOrigin=(25,25)
        self.SubGrids[1].RotateAngle=np.pi*3/4
        self.SubGrids[2].NewOrigin=(25,25)
        self.SubGrids[2].RotateAngle=np.pi/4
        self.SubGrids[3].NewOrigin=(25,25)
        self.SubGrids[3].RotateAngle=-np.pi/4
        
        return True
    
    def NeighborTest(self,GridID=0,Pts_end=(0.5,0)):
        '''Finding the neighbor subgrid and start point of a streamline
        '''
        #print('Original Grid:',GridID,Pts_end)
        #print(Pts_end[0],Pts_end[1])
        flag=False
        EdgeID=0
        if (Pts_end[0]>0.999999 and Pts_end[1]>0.0): #Streamline ended at out flow edge
            EdgeID=3
            flag=False
        elif (Pts_end[0]>0.0 and Pts_end[1]<0.000001): #Streamline ended at the bottom(b0) edge
            EdgeID=0
            flag=True
        elif (Pts_end[0]>0.0 and Pts_end[1]>0.999999): #Streamline ended at the top(b1) edge
            EdgeID=2
            flag=True
        
        NeighborGrid=0
        Neighbor_Init=(0,0)
        
        if (flag):#If streamline continue travel into other grid
            if (EdgeID==0): 
                NeighborGrid=self.NeighborID[GridID][0]
                Neighbor_Init=(Pts_end[0],1)
            if (EdgeID==2):
                NeighborGrid=self.NeighborID[GridID][1]
                Neighbor_Init=(Pts_end[0],0)
        
        #print('Neighbor Test:',flag,NeighborGrid,Neighbor_Init,EdgeID)
        
        return flag,NeighborGrid,Neighbor_Init,EdgeID
    
    def SLTrace(self,NSL=100,Pts=[]):
        """Trace Streamlines in the wellgrid using extended Pollock's Algorithm
           The streamline is traced from subgrids as sequence of 0,1,2,3
        
        Arguments
        ---------
        Pts - Initial particle start point at the well edge (alpha=0) (optional1)
        NSL - Number of Streamline (optional2)
        One has to made a choice of the two options

        Output
        ---------
        TOF        -- TOF array of of the [NSL] streamlines  
        SL         -- nodes array of of the [NSL] streamlines  
        TOF_end    -- a boundary TOF list of the [NSL] streamlines  
        SL_end     -- a boundary node list of the [NSL] streamlines 
        """
        TOF_end=[]
        SL_end=[]
        
        for i in range(4): #4 Subgrids
            
            if(len(Pts)==0):
                nsl=int(NSL*self.theta[i]/2/np.pi)
                Pts_init=PointOnUnitEdge(nsl) #Generating the start point along the well edge(alpha=0)
            else:
                nsl=len(Pts)
                Pts_init=Pts
                
            for j in range(nsl): #nsl streamlines
                GridID=i
                temp_trace=self.SubGrids[GridID].Trace1SL(Pts=Pts_init[j])
                
                SLtemp=RotateSL(temp_trace[3],Single=1,angle=self.SubGrids[GridID].RotateAngle)
                SLtemp=TranslateSL(SLtemp,Single=1,new_origin=self.SubGrids[GridID].NewOrigin)
                TOFtemp=temp_trace[5]
                
                flag=True
                while (flag==True): #the streamline will continue travel in another subgrid
                    Pts_end=temp_trace[2][-1]
                    temp_neighbor=self.NeighborTest(GridID,Pts_end) #test of crossing trace of a streamline
                    flag=temp_neighbor[0]
                    if(flag==True):
                        temp_trace=[]
                        SLtemp2=[]
                        TOFtemp2=[]
                        
                        GridID_next=temp_neighbor[1]
                        Pts_init_next=temp_neighbor[2]
                        #Pts and TOF base starts from previous node
                        temp_trace=self.SubGrids[GridID_next].Trace1SL(Pts=Pts_init_next,TOF_base=TOFtemp[-1])

                        SLtemp2=RotateSL(temp_trace[3],Single=1,angle=self.SubGrids[GridID_next].RotateAngle)
                        SLtemp2=TranslateSL(SLtemp2,Single=1,new_origin=self.SubGrids[GridID_next].NewOrigin)
                        TOFtemp2=temp_trace[5]
                        
                        #SLtemp=np.append(SLtemp,SLtemp2,axis=0)
                        #TOFtemp=np.append(TOFtemp,TOFtemp2,axis=0)
                        SLtemp=np.append(SLtemp,SLtemp2[1:],axis=0)
                        TOFtemp=np.append(TOFtemp,TOFtemp2[1:],axis=0)
                
                SL_end.append(SLtemp[-1])
                TOF_end.append(TOFtemp[-1])
                #Add all nodes and TOF into SL list
                self.SL.append(SLtemp)
                self.TOF.append(TOFtemp)
        
        
        #Plot the stremline
        plt.figure(figsize=(3, 3))
        plt.ylim(bottom=0,top=50)
        plt.xlim(left=0,right=50)
        plt.axes().set_aspect('equal')
        plt.title(r'Streamline in Physical Space ($x,y$)')
        
        #Grid edge
        Bound_vert=[self.Pts[0],self.Pts[1],self.Pts[2],self.Pts[3],self.Pts[0]]
        Internal_edge=[self.Pts[0],self.Pts[2],self.Pts[3],self.Pts[1]]
        
        plt.plot(*np.asarray(Bound_vert).T,lw=3,color='red')
        plt.plot(*np.asarray(Internal_edge).T,lw=2,ls='--',color='red')
        
        #Streamline
        for i in range(len(self.SL)):
            plt.plot(*np.asarray(self.SL[i]).T,lw=1,marker='o',markersize=0,color='blue')
        


        plt.show()
        return self.SL,self.TOF,SL_end,TOF_end