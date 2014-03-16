# RVIC Streamflow Routing Model

The RVIC Streamflow Routing model is a simple source to sink routing model.  The model represents each grid cell by a node in the channel network.  Unit hydrographs are developed that described the time distribution of flow from each source grid cell to a corresponding sink grid cell.  The development of the unit hydrographs is done as a pre process step (i.e. `make_parameters.py`).  The final step is the convolution of the unit hydrographs with fluxes from a land surface model, typically VIC (i.e. `rvic_model.py`).  

### Usage
See the [RVIC Wiki Page](https://github.com/jhamman/RVIC/wiki/RVIC-Wiki)

### Development History
1.  Based on the initial model of Lohmann, et al., 1996, Tellus, 48(A), 708-721
2.  Currently being coupled as a componenet model in the Comunity Earth System Model (RASM) as part of the Regional Arctic System Model Project.   

### Upcoming Development
1.  History file restarts to further support short model runs.
2.  Possible distributed application to support full grid routing (rather than source to sink only).
3.  Plans to couple directly to VIC.
