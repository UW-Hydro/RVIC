# Frequently Asked Questions

### Can the convolution be run on multiple processors?
At the moment, no.  However, many of the pieces needed to support this feature are already included in RVIC.  This is a feature that has been identified for a future version of RVIC.

### Can route flows to every grid cell in my domain using RVIC?
Yes, it is possible to route to every grid cell in the RVIC domain.  That said, it may not be the most efficient way to use RVIC.  RVIC is a source to sink routing model which means it doesn't track where streamflow is along its flow path.  When you route to every grid cell in the domain, you are duplicating a lot of calculations.  There are other routing models that route flow from grid cell to grid cell.

### Is there a C or Fortran Version of RVIC?
The RVIC convolution scheme has been added as a submodule extension to the VIC image driver, that submodule is written in C and allows RVIC to be coupled directly to VIC simulations. RVIC has also been coupled in the Community Earth System Model (CESM) as the streamflow routing model. For that project, a Fortran version of the "convolution" step was written.
