; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 10;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/MHD/Jet/gPLUTO_Comparison"
dir2 = "./"

title1 = "PLUTO"
title2 = "gPLUTO"

LOADCT,33
PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

imin=1
imax=4
; -- Load data from 1st directory --

PLOAD,/SIL, dir=dir1, nout
DISPLAY,/VBAR,x1=x1,x2=x2,/ysym,ALOG10(prs/rho),title=title1, $
        imin=imin,imax=imax

; -- Load data from 2nd directory --

PLOAD,/SIL, dir=dir2, nout
DISPLAY,/VBAR,x1=x1,x2=x3,/ysym,ALOG10(prs/rho),title=title2,nwin=1,$
        imin=imin,imax=imax

END

