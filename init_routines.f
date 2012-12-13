c     SUBROUTINES FOR INITIALIZATION
c     Read_NETCDF()

      SUBROUTINE READ_NETCDF( PATH, VAR , var_out)
      IMPLICIT NONE
      
      CHARACTER*50 PATH, VAR, var_out

      use netcdf
      implicit none
      integer                               :: ncId, VarId, &
                                               lonDimID, latDimId, timeDimId, &
                                               numLons, numLats, numTimes,    &
                                               status
      integer, dimension(nf90_max_var_dims) :: dimIDs
      real, dimension(:, :, :), allocatable :: var_out
      ...
      status = nf90_open(PATH, nf90_NoWrite, ncid)
      if(status /= nf90_NoErr) call handle_err(status)
      ...
      status = nf90_inq_varid(ncid, VAR, VarId)
      if(status /= nf90_NoErr) call handle_err(status)
c     What are the lengths of its constituent dimensions?
      status = nf90_inquire_variable(ncid, VarId, dimids = dimIDs)
      if(status /= nf90_NoErr) call handle_err(status)
      status = nf90_inquire_dimension(ncid, dimIDs(1), len = numLons)
      if(status /= nf90_NoErr) call handle_err(status)
      status = nf90_inquire_dimension(ncid, dimIDs(2), len = numLats)
      if(status /= nf90_NoErr) call handle_err(status)
      status = nf90_inquire_dimension(ncid, dimIDs(3), len = numTimes)
      if(status /= nf90_NoErr) call handle_err(status)
      allocate(var_out(numLons, numLats, numTimes))
      ...
      status = nf90_get_var(ncid, VarId, var_out)
      if(status /= nf90_NoErr) call handle_err(status)

      RETURN
      END

      SUBROUTINE WRITE_NETCDF (PATH, VAR)

      use netcdf
      implicit none
      
      CHARACTER*50 PATH, VAR, 
      REAL Qvals


      integer                               :: ncId, QVarId,status,          &
                                               timeDimId, &
                                               numTimes,    &
                                               i
      integer, dimension(nf90_max_var_dims) :: dimIDs
      ...
      status = nf90_open(PATH, nf90_Write, ncid)
      if(status /= nf90_NoErr) call handle_err(status)
      ...
      status = nf90_inq_varid(ncid, VAR, QVarId)
      if(status /= nf90_NoErr) call handle_err(status)
      ! How big is the netCDF variable, that is, what are the lengths of
      !   its constituent dimensions?
      status = nf90_inquire_variable(ncid, QVarId, dimids = dimIDs)
      if(status /= nf90_NoErr) call handle_err(status)
      status = nf90_inquire_dimension(ncid, dimIDs(3), len = numTimes)
      if(status /= nf90_NoErr) call handle_err(status)
      ...
      ! Make a temporary array the same shape as the netCDF variable.
      status = nf90_put_var(ncid, QVarId, Qvals)
      if(status /= nf90_NoErr) call handle_err(status)







