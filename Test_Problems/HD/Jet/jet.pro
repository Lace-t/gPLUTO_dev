; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 5;   Load this output
geom = "SPHERICAL"

LOADCT,/SIL,70
dir1 = "/Users/mignone/PLUTO/Test_Problems/HD/Jet/gPLUTO_Comparison"
dir2 = "/Users/mignone/gPLUTO-dev/Test_Problems/HD/Jet"

code_name1=(dir1.split('/'))[3]
code_name2=(dir2.split('/'))[3]

PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -- Load data from 1st directory --

PLOAD,/SIL, dir=dir1, nout
IF (geom EQ "SPHERICAL") THEN BEGIN
  pm = (Bx1^2 + Bx2^2 + Bx3^2)*0.5
  Q  = ALOG10(prs/pm)
  Q  = Bx3
  DISPLAY,x1=x1,x2=0.5*!PI-x2,/POLAR, /vbar, nwin=0, $
          title=code_name1, Q
ENDIF ELSE BEGIN
  DISPLAY,x1=x1,x2=x3,/vbar, vx2, nwin=0, title=code_name1,/ysym
ENDELSE

; -- Load data from 2nd directory --

PLOAD,/SIL, dir=dir2, nout
IF (geom EQ "SPHERICAL") THEN BEGIN
  pm = (Bx1^2 + Bx2^2 + Bx3^2)*0.5
  Q  = Bx3; ALOG10(prs/pm)
  DISPLAY,x1=x1,x2=0.5*!PI-x2,/POLAR, /vbar, nwin=2, $
          title=code_name2, Q
ENDIF ELSE BEGIN
  DISPLAY,x1=x1,x2=x3,/vbar, vx2, nwin=1, title=code_name2,/ysym

ENDELSE

END

