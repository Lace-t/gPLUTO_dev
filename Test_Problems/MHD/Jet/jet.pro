PLOAD, /SIL, 0

LOADCT,/SIL,70
REVERSE_COLORMAP
FOR nout = 0,nlast DO BEGIN
  PLOAD, /SIL, nout
  DISPLAY, /ysym, /VBAR, x1=x1, x2=x2, alog10(rho), $
           title = "t = "+ARG2STR(t[nout]), nwin=0
  WAIT,0.25
ENDFOR


END
