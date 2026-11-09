dir1 = "Serial"
dir2 = "nprocs_2"
PLOAD,dir=dir1,/NODATA,0

FOR nout = 0, nlast DO BEGIN

  PLOAD,dir=dir1,/SIL, nout
  rho1 = rho

  PLOAD,dir=dir2,/SIL, nout
  rho2 = rho

  PRINT,"nout = ",nout,"; MAX(ABS(rho1-rho2)) = ", MAX(ABS(rho1-rho2))
ENDFOR

END


