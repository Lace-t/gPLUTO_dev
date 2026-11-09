; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 100;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/MHD/Torus/gPLUTO_Comparison"
dir2 = "/Users/mignone/gPLUTO-dev/Test_Problems/MHD/Torus"

title1 = "PLUTO"
title2 = "gPLUTO"

code_name1=(dir1.split('/'))[3]
code_name2=(dir2.split('/'))[3]

PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

LOADCT,/SIL,33
; -----------------------------------------------
; 1. Load data from 1st directory 
; -----------------------------------------------

PLOAD,/SIL, dir=dir1, nout 
;DISPLAY, x1=x1,x2=0.5*!PI - x2, /POLAR, /VBAR, ims=[400,800], title=title1,bx1
DISPLAY, x1=x1,x2=x3, /VBAR, ims=[400,800], title=title1,bx2

; -----------------------------------------------
; 2. Load data from 2nd directory 
; -----------------------------------------------

PLOAD,/SIL, dir=dir2, nout
;DISPLAY, x1=x1,x2=0.5*!PI- x2, /POLAR, /VBAR, ims=[400,800], title=title2, $
;         nwin=1,bx1
DISPLAY, x1=x1,x2=x3,/VBAR, ims=[400,800], title=title2, $
         nwin=1,bx2


END

