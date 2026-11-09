; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 10;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/HD/Stellar_Wind/gPLUTO_Comparison"
dir2 = "./"

title1 = "PLUTO   <---> gPLUTO"

LOADCT,13
PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -- Load data from 1st directory --
PLOAD,/SIL, dir=dir1, /VTK, nout

IF (nx2 EQ 1) THEN BEGIN
  PRINT,'2D Data'
  PLOAD,/SIL, dir=dir1, /VTK, nout
  Q = FLTARR(2*nx1,nx3)

  Q[nx1-1:0:-1,*] = alog10(REFORM(rho))

  ; -- Load data from 2nd directory --
  PLOAD,/SIL,dir=dir2, /VTK, nout
  Q[nx1:2*nx1-1,*] = alog10(REFORM(rho))

  DISPLAY, /vbar, Q, title=title1
ENDIF

IF (nx2 GT 1) THEN BEGIN    ; 3D Cart
  PRINT,'3D Data'
  PLOAD,/SIL, dir=dir1, /VTK, nout
  DISPLAY,/vbar,alog10(rho[*,*,nx3/2]),title="PLUTO",ims=4

  ; -- Load data from 2nd directory --
  PLOAD,/SIL,dir=dir2, /VTK, nout
  
  DISPLAY,/vbar,alog10(rho[*,*,nx3/2]),nwin=1,title="gPLUTO",ims=4
ENDIF

END

