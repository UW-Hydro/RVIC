# RVIC Streamflow Routing Model

The RVIC Streamflow Routing model is a simple source to sink routing model.  The model represents each grid cell by a node in the channel network.  Unit hydrographs are developed that described the time distribution of flow from each source grid cell to a corresponding sink grid cell.  The development of the unit hydrographs is done as a pre process step (i.e. `make_parameters.py`).  The final step is the convolution of the unit hydrographs with fluxes from a land surface model, typically VIC (i.e. `rvic_model.py`).  

### Usage
Examples of all configuration files can be found in RVIC/config/ directory

#### Make domain file from arcinfo fraction file

`./fraction2domain.bash fraction_file.asc`

#### Make rvic parameter file from existing uhs files 
*** note that only C program format files are currently supported

`./uhs2paramfile.py config_file.cfg`

#### Make rvic parameter file from pour points

`./make_parameters.py config_file.cfg -np $N_PROCESSORS`

where `$N_PROCESSORS` is the number of processors to be used.  This uses the Python [multiprocessing](http://docs.python.org/2/library/multiprocessing.html) module.  

#### Run the rvic model convolution
`./rvic_model.py config_file.cfg`

### Development History
1.  Based on the initial model of Lohmann, et al., 1996, Tellus, 48(A), 708-721
2.  Currently being coupled as a componenet model in the Comunity Earth System Model (RASM) as part of the Regional Arctic System Model Project.   

### Upcoming Development
1.  Write core convolution routine in C for speed.  
2.  History file restarts to further support short model runs.
3   Possible distributed application to support full grid routing (rather than source to sink only).  
3.  Plans to couple directly to VIC
