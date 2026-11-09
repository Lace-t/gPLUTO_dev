; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 1;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/RMHD/Rotor/gPLUTO_Comparison"
dir2 = "/Users/mignone/gPLUTO-dev/Test_Problems/RMHD/Rotor"

code_name1=(dir1.split('/'))[3]
code_name2=(dir2.split('/'))[3]

PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -----------------------------------------------
; 1. Load data from 1st directory 
; -----------------------------------------------

LOADCT,/SIL,13
PLOAD,/SIL, dir=dir1, nout 
DISPLAY, x1=x1,x2=x2, /POLAR, /VBAR, ims=[400,400], title=dir1, alog10(prs)

; -----------------------------------------------
; 2. Load data from 2nd directory 
; -----------------------------------------------

PLOAD,/SIL, dir=dir2, nout
DISPLAY, x1=x1,x2=x2, /POLAR, /VBAR, ims=[400,400], title=dir2, nwin=1, alog10(prs)

END

