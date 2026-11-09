; ----------------------------------------------------------------------
; Visual comparison for the jet propagation problem
;
; Last Modified:  26 Fen 2025 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 1;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/RMHD/Jet/gPLUTO_Comparison"
dir2 = "/Users/mignone/gPLUTO-dev/Test_Problems/RMHD/Jet"

title1 = "PLUTO"
title2 = "gPLUTO"

code_name1=(dir1.split('/'))[3]
code_name2=(dir2.split('/'))[3]

PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -----------------------------------------------
; 1. Load data from 1st directory 
; -----------------------------------------------

LOADCT,/SIL,13
PLOAD,/SIL, dir=dir1, nout 
DISPLAY, x1=x1,x2=x3, /VBAR, title=title1, alog10(rho)

; -----------------------------------------------
; 2. Load data from 2nd directory 
; -----------------------------------------------

PLOAD,/SIL, dir=dir2, nout
DISPLAY, x1=x1,x2=x3, /VBAR, title=title2, nwin=1, alog10(rho)

END

