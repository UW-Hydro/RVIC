# Frequently Asked Questions

### Can the convolution be run on multiple processors?
At the moment, no.  However, many of the pieces needed to support this feature are already included in RVIC.  This is a feature that has been identified for a future version of RVIC.

### Can route flows to every grid cell in my domain using RVIC?
Yes, it is possible to route to every grid cell in the RVIC domain.  That said, it may not be the most efficient way to use RVIC.  RVIC is a source to sink routing model which means it doesn't track where streamflow is along its flow path.  When you route to every grid cell in the domain, you are duplicating a lot of calculations.  There are other routing models that route flow from grid cell to grid cell.

### Is there a C or Fortran Version of RVIC?
RVIC has been coupled in the Community Earth System Model (CESM) as the streamflow routing model.  For that project, a Fortran version of the "convolution" step was written.  Work is currently underway to couple RVIC to the stand-alone image driver in VIC 5.0 ([VIC231](https://github.com/UW-Hydro/VIC/pull/231)).
