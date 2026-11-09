
dir1 = "Ref"
dir2 = "."
FOR nout = 0,10 DO BEGIN
  PLOAD,/SIL,dir=dir1,nout
  rho1 = rho

  PLOAD,/SIL,dir=dir2,nout
  rho2 = rho

  PRINT,"nout = ",nout,"; MAX(drho) = ", MAX(ABS(rho1-rho2))
ENDFOR


END

