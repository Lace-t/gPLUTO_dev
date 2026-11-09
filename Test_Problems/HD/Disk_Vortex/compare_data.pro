; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 1;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/HD/Disk_Vortex/gPLUTO_Comparison"
dir2 = "./"

LOADCT,13
PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -- Load data from 1st directory --
PLOAD,/SIL, dir=dir1, nout
DISPLAY,x1=x1,x2=x2,/POLAR, rho, title="PLUTO",ims=[400,400],/VBAR


; -- Load data from 2nd directory --
PLOAD,/SIL,dir=dir2, nout
DISPLAY,x1=x1,x2=x2,/POLAR, rho, title="gPLUTO-dev", nwin=1,ims=[400,400],/VBAR

END

