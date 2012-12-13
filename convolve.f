c     Beta Fortran convolution routine for implementation within RASM
c     Step 0 - Setup
c     Step 1 - Load IRF NetCDFs (Fraction and IRF)
c     Step 2 - Load Runoff/Baseflow Fluxes
c     Step 3 - Load Grid Cell Areas
c     Step 4 - Make Convolution

      PROGRAM MAKE_CONVOLUTION

c     Step 0 - Setup
      
      IMPLICIT NONE

      CHARACTER*50 RCSID

      INTEGER     N, I, J, DAYS, NDAY, II, JJ  
      INTEGER     NCOL,NROW,ICOL,NOB,PMAX,KE,UH_DAY
      INTEGER     CATCHIJ(PMAX,2)
      INTEGER     NYR
      REAL        UH_S(PMAX,KE+UH_DAY-1)
      REAL        BASE(DAYS), RUNO(DAYS), FLOW(DAYS) 
      REAL        FRACTION(NCOL,NROW)

      REAL        PI, RERD, FACTOR, FACTOR_SUM

      LOGICAL TORF

      PARAMETER   (RERD  = 6371229.0)    !radius of earth in meters

      CHARACTER*20 LOC
      CHARACTER*72 INPATH

      INTEGER DPREC, CLEN

      REAL        JLOC, ILOC
      REAL        XC, YC, SIZE
      REAL        AREA, AREA_SUM

      INTEGER     IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER MO(12*NYR),YR(12*NYR)

      PI = ATAN(1.0) * 4.0
      
      use netcdf
 
c     Step 1 - Load IRF NetCDF (Fraction and IRF)
      status = nf90_open(path = "/raid/jhamman/temp_uh_files/run8/DVIUP_RASM_UH.nc",
     & cmode = nf90_nowrite, ncid = fraction)
      if (status /= nf90_noerr) call handle_err(status)
   

c     Step 2 - Load Runoff/Baseflow Fluxes


c     Step 3 - Load Grid Cell Areas


c     Step 4 - Make Convolution
