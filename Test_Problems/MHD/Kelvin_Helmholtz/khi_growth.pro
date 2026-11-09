
; -- Evolution plot --
PLOAD,/SIL,0

LOADCT,/SIL,13
FOR nout = 0,-nlast DO BEGIN
  PLOAD,/SIL,nout
  DISPLAY, x1=x1, x2=x2, /vbar, tr1, ims=[200,400],$
           title = "t = "+ARG2STR(t[nout])
  WAIT,0.2
ENDFOR

; -- Growth rate --
LOADCT,/SIL, 0
WINDOW,1
dv = FLTARR(nlast+1)

Q = READ_ASCII("khi.dat",comm="#")
Q = Q.(0)
t = Q[0,*]
dv = Q[1,*]

; *************************************
; Scale exact growth rate
;
; t0 = a/cs
; (ka)*(x/a) - (w*t0)*(t/t0)
; w*a/cs = wc --> w = wc*cs/a
; *************************************

a  = 0.02
cs = 1.0
om = 0.071*cs/a

; -- Plot numerical and analytical growth rate --
PLOT,  t, dv, /YLOG, psym=-2, position=[0.15,0.15,0.95,0.95],$
       xtitle="t",ytitle = "Growth Rate", charsize=2., back=255,col=0
OPLOT, t, 1.e-6*exp(om*t),line=2,col=0

LEGEND2, ["e!U!4x!3t!N","!4D!3v!Dy!N(t)"], line=[2,0],psym=[0,-2],$
         col=0,textcol=0,/BOT,/RIGHT, charsize=2.5, charthick=2

END

