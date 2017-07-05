
Source Code of Nearwellbore Streamline Tracing (SPE-182614-MS)- v0.2.0
==============================================================================================

Bin Wang (binwang.0213@gmail.com), Yin Feng

Department of Petroleum Engineering, Univeristy of Louisiana at Lafayette, Lafayette, US, 70506

Corresponding author: Yin Feng, yin.feng@louisiana.edu
![Image of Embedded Method](https://github.com/BinWang0213/Nearwellbore_Streamline/blob/master/images/Embedded_Field.png)

This code includes following three algorithm to construct the streamlines within a 2.5D wellblock
1. **Standard Subdivision Method** - The wellblock is subdivided by four triangular sub-grids and streamlines are traced using extended Pollock's algorithm ( Chapter 8 on Datta-Gupta and King, 2007)
2. **Embedded Method** - The velocity field is solved using virtual-boundary-element method (VBEM) and streamlines are traced using numercial integertor.The detailed description can be found on SPE-182614-MS
3. **Simple Fill-Grid Method** - The average time-of-flight (TOF) is estimated and no streamlines are explicited traced (Natvig et al. 2007). 

Here is a nicely formatted [ipython notebook version](https://github.com/BinWang0213/Nearwellbore_Streamline/blob/master/Example-QuickStart.ipynb) of running the algorithm with Case 3 in SPE-182614-MS. 

Reference:  
[1] Datta-Gupta, A. and King, M.J., 2007. Streamline simulation: theory and practice. Textbook Series, SPE, Richardson, Texas 11.  
[2] Natvig, J.R., Lie, K.A., Eikemo, B. and Berre, I., 2007. An efficient discontinuous Galerkin method for advective transport in porous media. Advances in water resources, 30(12), 2424-2438.  doi:10.1016/j.advwatres.2007.05.015  

Citation
--------

Wang, B., Feng, Y., Du, J., et al. (2017) An Embedded Grid-Free Approach for Near Wellbore Streamline Simulation. doi:10.2118/SPE-182614-MS

License
-------
This code is released under the terms of the BSD license, and thus free for commercial and research use. Feel free to use the code into your own project with a PROPER REFERENCE.  
