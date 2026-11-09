; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 8;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/HD/Disk_Planet/gPLUTO_Comparison"
dir2 = "/Users/mignone/gPLUTO-dev/Test_Problems/HD/Disk_Planet"

code_name1=(dir1.split('/'))[3]
code_name2=(dir2.split('/'))[3]

PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -----------------------------------------------
; 1. Load data from 1st directory 
; -----------------------------------------------

LOADCT,/SIL,13
PLOAD,/SIL, dir=dir1, /VTK, nout 
Q = alog(rho[*,nx2-1,*])
DISPLAY, x1=x1,x2=x3, /POLAR, /VBAR, ims=[400,400], title=dir1, Q

; -----------------------------------------------
; 2. Load data from 2nd directory 
; -----------------------------------------------

PLOAD,/SIL, dir=dir2, /VTK, nout
Q = alog(rho[*,nx2-1,*])
DISPLAY, x1=x1,x2=x3, /POLAR, /VBAR, ims=[400,400], title=dir2, nwin=1, Q

END

