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
from .geometry import *

                            # Fill-Grid Method #

def Fillgrid(Pts_e=[],Pts_w=[],Qw=[],h=26.25,phi=0.2):
    """Calculate the average time-of-flight(TOF) based on the fill-grid method
        
        Arguments
        ---------
        Qw     -- Wellbore flow rate, stb/day
        Pts_e  -- Vertices of boundary edge e.g. [(0,0),(0,1),(1,1),(1,0)]
        h      -- Wellblock thickness
        phi    -- Wellblock porosity
    
        Bulk_Vol -- Rock bulk volume
        Qw_total -- Total wellbore injection/production rate
        TOF_Avg  -- Average TOF for the wellblock
    """
    Bulk_Vol=PolygonArea(Pts_e)*h
    Qw_total=sum(Qw)

    TOF_Avg=Bulk_Vol*phi/Qw_total/5.6145
    #Plot
    domain_min=(min(np.asarray(Pts_e)[:,0]),min(np.asarray(Pts_e)[:,1]))
    domain_max=(max(np.asarray(Pts_e)[:,0]),max(np.asarray(Pts_e)[:,1]))
    
    plt.figure(figsize=(3, 3))
    plt.axes().set(xlim=[domain_min[0], domain_max[0]],
                   ylim=[domain_min[1], domain_max[1]],aspect='equal')

    #Domain boundary
    plt.plot(*np.asarray(Pts_e).T,lw=1,color='black')
    plt.scatter(*np.asarray(Pts_w).T,s=20,color='red')
    
    plt.axes().set_title('Fill-Grid Method\n $TOF_{avg}=%.4f day$'%(TOF_Avg))
    plt.show()
    return TOF_Avg


