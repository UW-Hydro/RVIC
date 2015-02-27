# Frequently Asked Questions

### Can the convolution be run on multiple processors?
At the moment, no.  However, many of the pieces needed to support this feature are already included in RVIC.  This is a feature that has been identified for a future version of RVIC.

### Can route flows to every grid cell in my domain using RVIC?
Yes, it is possible to route to every grid cell in the RVIC domain.  That said, it may not be the most efficent way to use RVIC.  RVIC is a source to sink routing model which means it doesn't track where streamflow is along its flow path.  When you route to every grid cell in the domain, you are duplicating a lot of calculations.  There are other routing models that route flow from grid cell to grid cell.

### Is there a C or Fortran Version of RVIC?
RVIC has been coupled in the Community Earth System Model (CESM) as the streamflow routing model.  For that project, a Fortran version of the "convolution" step was written.  At this time, there is not a C version of this routing model.  In the future, a C binding may be created to coupled with the stand-alone VIC model version 5.
