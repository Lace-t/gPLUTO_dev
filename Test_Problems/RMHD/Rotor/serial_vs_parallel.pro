; ----------------------------------------------------------------------
; Check differences between data produced with serial_vs_parallel.sh
; Last Modified: 17 Jun 2024 by A. Mignone (andrea.mignone@unito.it)
; ----------------------------------------------------------------------
nout = 1;   Load this output

dir1 = "Serial"
dir2 = "nprocs_8"

PRINT,"dir1 = ",dir1
PRINT,"dir2 = ",dir2

; -- Load data from 1st directory --

PLOAD,/SIL, dir=dir1, nout
FOR nout = 0,nlast DO BEGIN
  PLOAD,/SIL, dir=dir1, nout
  rho_1 = rho
  Bx1_1 = Bx1
  Bx2_1 = Bx2
  Bx3_1 = Bx3
  IF (ISA(Bx1s)) THEN Bx1s_1 = Bx1s
  IF (ISA(Bx2s)) THEN Bx2s_1 = Bx2s
  IF (ISA(Bx3s)) THEN Bx3s_1 = Bx3s
  IF (ISA(Ex1s)) THEN Ex1s_1 = Ex1s
  IF (ISA(Ex2s)) THEN Ex2s_1 = Ex2s
  IF (ISA(Ex3s)) THEN Ex3s_1 = Ex3s

  ; -- Load data from 2nd directory --
  PLOAD,/SIL, dir=dir2, nout
  rho_2 = rho
  Bx1_2 = Bx1
  Bx2_2 = Bx2
  Bx3_2 = Bx3
  IF (ISA(Bx1s)) THEN Bx1s_2 = Bx1s
  IF (ISA(Bx2s)) THEN Bx2s_2 = Bx2s
  IF (ISA(Bx3s)) THEN Bx3s_2 = Bx3s
  IF (ISA(Ex1s)) THEN Ex1s_2 = Ex1s
  IF (ISA(Ex2s)) THEN Ex2s_2 = Ex2s
  IF (ISA(Ex3s)) THEN Ex3s_2 = Ex3s

  ; -- Print data difference --

  PRINT,"drho = ", MAX(ABS(rho_1 - rho_2))
  PRINT,"dBx  = ", MAX(ABS(Bx1_1 - Bx1_2))
  PRINT,"dBy  = ", MAX(ABS(Bx2_1 - Bx2_2))
  PRINT,"dBz  = ", MAX(ABS(Bx3_1 - Bx3_2))
  IF (ISA(Bx1s)) THEN PRINT,"dBxs  = ", MAX(ABS(Bx1s_1 - Bx1s_2))
  IF (ISA(Bx2s)) THEN PRINT,"dBys  = ", MAX(ABS(Bx2s_1 - Bx2s_2))
  IF (ISA(Bx3s)) THEN PRINT,"dBzs  = ", MAX(ABS(Bx3s_1 - Bx3s_2))

  IF (ISA(Ex1s)) THEN PRINT,"dExs  = ", MAX(ABS(Ex1s_1 - Ex1s_2))
  IF (ISA(Ex2s)) THEN PRINT,"dEys  = ", MAX(ABS(Ex2s_1 - Ex2s_2))
  IF (ISA(Ex3s)) THEN PRINT,"dEzs  = ", MAX(ABS(Ex3s_1 - Ex3s_2))

  ; -- PLOT --

  DISPLAY,x1=x1,x2=x2,/POLAR,ABS(rho_1 - rho_2) GT 1.e-15, title = "t = "+ARG2STR(nout)
  WAIT,0.5
ENDFOR

END

