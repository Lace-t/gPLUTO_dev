; ----------------------------------------------------------------------
; IDL script to visualize data for the MHD blast problem.
;
; Last Modified: Oct 22, 2018 by A. Mignone
; ----------------------------------------------------------------------

FUNCTION PlutoVar, name  ; Set the variable to be plotted

  COMMON PLUTO_VAR

;  name = "!4q!3"  & RETURN, rho
  name = "p"      & RETURN, prs
END

; -----------------------------------------------
; 0. Define directory,
;    enable/disable image writing
; -----------------------------------------------
dir        = "./"
save_image = 1
fbase_name = "bmag0.01"

PLOAD, dir=dir,0
B0 = MAX(bx1)
!X.STYLE=1
!Y.STYLE=1
!P.CHARSIZE=1.8
!P.CHARTHICK=2

; -----------------------------------------------
; 1. Get max and min values of the quantity
;    to be plotted.
;    The plasma beta is computed from the
;    initial condition.
; -----------------------------------------------
maxq = -1.e3
minq = -maxq
FOR nout = 0, nlast DO BEGIN
  PLOAD, dir=dir,nout,/SILENT
  beta = 2.0*prs/(bx1^2 + bx2^2)
  q    = PlutoVar(vname)
  IF (nout EQ 0) THEN BEGIN
    beta_title = "!4b!3 = "+ARG2STR(MIN(beta),FORMAT='(F6.2)')
  ENDIF

  maxq = MAX([MAX(q),maxq])
  minq = MIN([MIN(q),minq])
ENDFOR

; -----------------------------------------------
; 2. Loop over datafiles and plot
; -----------------------------------------------
FOR nout = 0, nlast DO BEGIN
  PLOAD, dir=dir,nout,/SILENT
  LOADCT,33,/SIL
  beta = 2.0*prs/(bx1^2 + bx2^2)
  q = PlutoVar(vname) 
  title = "MHD Blast ("+ARG2STR(beta_title)+"), "+vname+"(t = "+ARG2STR(t[nout])+")
  DISPLAY,x1=x1,x2=x2, q, /VBAR, title=title,label1="x",label2="y",$
          imax = maxq*0.9, imin=minq*1.1, ims=[512,512]
     
  nlev   = 20
  levels = MIN(Ax3) + findgen(nlev)*(MAX(Ax3) - MIN(Ax3))/(nlev-1.0)
  LOADCT,0,/SIL
  CONTOUR,ax3,x1,x2,/overplot,/closed,levels=levels,col=240
  WAIT,0.1
 
  IF (save_image) THEN BEGIN
    fname = fbase_name+"."+STRCOMPRESS(STRING(nout,format='(I4.4)'),/REMOVE_ALL)
    GET_FRAME,name=fname,/gif
  ENDIF
  
ENDFOR

END

