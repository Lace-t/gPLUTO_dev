dir = ["run_16","run_32","run_64", "run_128", "run_256"]
err = FLTARR(16)

ndir = (size(dir))[1] 
FOR i = 0, ndir-1 DO BEGIN

  ; -- Load initial sol --
  PLOAD,dir=dir[i], /SIL,0
  dV   = dx1#dx2
  rho0 = rho

  ; -- Load final sol --
  PLOAD,dir=dir[i], /SIL,1
  dV   = dx1#dx2

  ; -- Compute error --
  err[i] = TOTAL(dV*ABS(rho-rho0))/TOTAL(dV)

  order_str = "-" 
  IF (i GT 0) THEN BEGIN  
    order     = alog10(err[i]/err[i-1])/alog10(0.5)
    order_str = STRING(order)
  ENDIF
  PRINT, "Dir = "+dir[i]+"; err = ",err[i],"; order = "+order_str

ENDFOR


END
