
nout = 2
dir1 = "./Full_Run"
dir2 = "./"

; -------------------------------------
; Load data from 1st directory
; -------------------------------------

PLOAD,dir=dir1,/SIL,nout
rho_01  = rho
prs_01  = prs
IF (ISA(Bx1s)) THEN Bx1s_01 = Bx1s 
IF (ISA(Bx2s)) THEN Bx2s_01 = Bx2s 
IF (ISA(flag)) THEN flag_01 = flag

; -------------------------------------
; Load data from 2nd directory
; -------------------------------------

PLOAD,dir=dir2,/SIL,nout
rho_02  = rho
prs_02  = prs
IF (ISA(Bx1s)) THEN Bx1s_02 = Bx1s 
IF (ISA(Bx2s)) THEN Bx2s_02 = Bx2s 
IF (ISA(flag)) THEN flag_02 = flag

; -------------------------------------
; Compute difference
; -------------------------------------

PRINT,"err(rho)",  MAX(ABS(rho_01-rho_02))
PRINT,"err(prs)",  MAX(ABS(prs_01-prs_02))
IF (ISA(Bx1s)) THEN PRINT,"err(Bx1s)", MAX(ABS(Bx1s_01-Bx1s_02))
IF (ISA(Bx2s)) THEN PRINT,"err(Bx2s)", MAX(ABS(Bx2s_01-Bx2s_02))
IF (ISA(flag)) THEN PRINT,"err(flag)", MAX(ABS(flag_01-flag_02))

END
