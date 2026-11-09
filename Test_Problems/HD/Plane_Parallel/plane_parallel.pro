; ----------------------------------------------------------------------
; Check binary equality of two different data set
;
; Last Modified:  02 Jul 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 1;   Load this output

dir1 = "/Users/mignone/PLUTO/Test_Problems/HD/Spherical/Plane_Parallel/gPLUTO_Comparison"
dir2 = "./"

LOADCT,13
PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -- Load data from 1st directory --
PLOAD,/SIL, dir=dir1, 0
rho_0 = rho
vx1_0 = vx1
vx2_0 = vx2
prs_0 = prs

PLOAD,/SIL, dir=dir1, nout
DISPLAY, x1=x1,x2=x2, /vbar, prs, title=dir1

PRINT,"err rho (dir = "+dir1+") = ", TOTAL(ABS(rho-rho_0))/nx1/nx2
PRINT,"err vx1 (dir = "+dir1+") = ", TOTAL(ABS(vx1-vx1_0))/nx1/nx2
PRINT,"err vx2 (dir = "+dir1+") = ", TOTAL(ABS(vx2-vx2_0))/nx1/nx2
PRINT,"err prs (dir = "+dir1+") = ", TOTAL(ABS(prs-prs_0))/nx1/nx2

; -- Load data from 2nd directory --
PLOAD,/SIL, dir=dir2, 0
rho_0 = rho
vx1_0 = vx1
vx2_0 = vx2
prs_0 = prs

PLOAD,/SIL, dir=dir2, nout

DISPLAY, x1=x1,x2=x2,/vbar, prs,nwin=1, title=dir2

PRINT,"err rho (dir = "+dir2+") = ", TOTAL(ABS(rho - rho_0))/nx1/nx2
PRINT,"err vx1 (dir = "+dir2+") = ", TOTAL(ABS(vx1 - vx1_0))/nx1/nx2
PRINT,"err vx2 (dir = "+dir2+") = ", TOTAL(ABS(vx2 - vx2_0))/nx1/nx2
PRINT,"err prs (dir = "+dir2+") = ", TOTAL(ABS(prs - prs_0))/nx1/nx2


END

