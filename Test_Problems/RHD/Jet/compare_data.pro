; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------

nout = 6;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/RHD/Jet/gPLUTO_Comparison"
dir2 = "./"

title1 = "PLUTO   <---> gPLUTO"

LOADCT,33
PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -- Load data from 1st directory --

PLOAD,/SIL, dir=dir1, /vtk, nout
Q = FLTARR(2*nx1,nx2)

Q[nx1-1:0:-1,*] = ALOG10(REFORM(rho))

; -- Load data from 2nd directory --
PLOAD,/SIL,dir=dir2, /vtk, nout
Q[nx1:2*nx1-1,*] = alog10(REFORM(rho))


DISPLAY, /vbar, Q, title=title1

WAIT,0.25
;DISPLAY, x1=x1,x2=x3,/vbar, alog10(prs),nwin=1, title=title2

END

