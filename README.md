# RVIC Streamflow Routing Model

[![Build Status](https://travis-ci.org/jhamman/RVIC.svg?branch=master)](https://travis-ci.org/jhamman/RVIC)

The RVIC Streamflow Routing model is a simple source to sink routing model.  The model represents each grid cell by a node in the channel network.  Unit hydrographs are developed that described the time distribution of flow from each source grid cell to a corresponding sink grid cell.  The development of the unit hydrographs is done as a pre process step (i.e. `rvic parameters`).  The final step is the convolution of the unit hydrographs with fluxes from a land surface model, typically VIC (i.e. `rvic convolution`).  

### Usage
See the [RVIC Wiki Page](https://github.com/jhamman/RVIC/wiki/RVIC-Wiki)

### Development:
1.  Based on the initial model of Lohmann, et al., 1996, Tellus, 48(A), 708-721
2.  Currently being coupled as a componenet model in the Comunity Earth System Model (RASM) as part of the Regional Arctic System Model Project.

More information on the development and testing of the RVIC model can be found [here](https://github.com/jhamman/RVIC/wiki/Development-and-Testing).

### References:
**RVIC Primary Reference**

1.  Hamman, J.J., and Coauthors.  (Manuscript in preparation - 2014).  A simple coupled streamflow routing model for global and regional climate models.

**Primary Historical Reference**

2.  Lohmann, D., R. Nolte-Holube, and E. Raschke, 1996: A large-scale horizontal routing model to be coupled to land surface parametrization schemes, Tellus, 48(A), 708-721.

3.  Lohmann, D., E. Raschke, B. Nijssen and D. P. Lettenmaier, 1998: Regional scale hydrology: I. Formulation of the VIC-2L model coupled to a routing model, Hydrol. Sci. J., 43(1), 131-141.

**Other Historical References**

4.  Nijssen, B.N., D.P. Lettenmaier, X. Liang, S.W. Wetzel, and E.F. Wood, 1997: Streamflow simulation for continental-scale river basins, Water Resour. Res., 33(4), 711-724.

5.  Lohmann, D., E. Raschke, B. Nijssen, and D. P. Lettenmaier, 1998: Regional Scale Hydrology: II. Application of the VIC-2L Model to the Weser River, Germany, Hydrological Sciences Journal, 43(1), 143-157.
