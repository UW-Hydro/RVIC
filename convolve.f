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
      CALL READ_NETCDF("/raid/jhamman/temp_uh_files/run8/DVIUP_RASM_UH.nc","fraction","FRAC")
      CALL READ_NETCDF("/raid/jhamman/temp_uh_files/run8/DVIUP_RASM_UH.nc","unit_hydrograph", "IRF")


c     Step 2 - Load Runoff/Baseflow Fluxes
      CALL READ_NETCDF("/nfs/thermal/raid/jhamman/RASM_results/r33/r33RBVIC70.vic.ha.1990s.qflux.nc","Runoff","Runoff")
      CALL READ_NETCDF("/nfs/thermal/raid/jhamman/RASM_results/r33/r33RBVIC70.vic.ha.1990s.qflux.nc","Baseflow","Baseflow")
      FLUX = Rufoff+Baseflow


c     Step 3 - Load Grid Cell Areas
      CALL READ_NETCDF("/nfs/hydro6/raid/nijssen/rasm/masks/racm_masks_121108.nc","DOMA_AREA","AREA")
      AREA = AREA*RERD*RERD


c     Step 4 - Make Convolution
      n = len(FLUX)
      m = len(IRF)
      do i=1,n
         do j=1,m
            T_STEP = i+j
            Q(T_STEP) = Q(T_STEP)+SUM(FLUX(i,:,:)*IRF(J,:,:)*AREA(:,:)*FRAC(:,:),mask=FRAC.gt.0)
         end do
      end do


c     Step 5 - Write to netCDF
      CALL WRITE_NETCDF("/out.nc",Q)
      
