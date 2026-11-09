nout = 0
dir1 = "Ref"
dir2 = "./"


LOADCT,13
PLOAD,/SIL,dir=dir1, nout
Bx1s_ref = Bx1s
Bx2s_ref = Bx2s
DISPLAY,Bx2s_ref,/VBAR,ims=.5

PLOAD,/SIL,dir=dir2, nout
DISPLAY,Bx2s,/VBAR,nwin=1,ims=.5


PRINT,MAX(ABS(Bx1s_ref - Bx1s))
PRINT,MAX(ABS(Bx2s_ref - Bx2s))

;PLOT,Bx1s_ref[*,nx2/2] 
;OPLOT,Bx1s[*,nx2/2] ,psym=1

END

