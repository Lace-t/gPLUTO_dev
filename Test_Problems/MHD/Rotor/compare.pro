nout = 0
PLOAD,/SIL,dir="Ref",nout
Bx1_r = Bx1s & Bx2_r = Bx2s

PLOAD,/SIL,dir=".",nout

PLOT,  Bx1_r[nx2/2,*]
OPLOT, Bx1s [nx2/2,*],psym=1

PRINT,"err = ",MAX(ABS(Bx1_r-BX1s))
; OPLOT,x1, Bx2_r, color=200
; OPLOT,x1, Bx2, psym=1, color=200

END


