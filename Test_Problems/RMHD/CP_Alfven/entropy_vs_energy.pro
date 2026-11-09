nout = 10
dir1 = "Energy"
dir2 = "Entropy"
PLOAD, dir=dir1, /SIL, 0

FOR nout = 0, nlast DO BEGIN

  ; Read energy data
  PLOAD, dir=dir1, /SIL, nout
  Bx_1  = Bx1 
  By_1  = Bx2
  Bz_1  = Bx3
  prs_1 = prs
  pm_1  = Bx1^2 + Bx2^2 + Bx3^2
  Em_1  = Ex1^2 + Ex2^2 + Ex3^2

  ; Read entropy data
  PLOAD, dir=dir2, /SIL, nout
  Bx_2  = Bx1
  By_2  = Bx2
  Bz_2  = Bx3
  prs_2 = prs
  pm_1  = Bx1^2 + Bx2^2 + Bx3^2
  Em_2  = Ex1^2 + Ex2^2 + Ex3^2

  ; Plot (mag. field)
  strtime = "time = "+ARG2STR(t[nout])
  PLOT, x1, By_1, charsize=1.4, $
        title = strtime, xtitle="x", ytitle="By"; , yrange=[1.4,2]
  OPLOT,x1, By_2, psym=1

  ; PLOT, x1, ABS(prs_1/prs_2-1.0), charsize=1.4, $
  ;       title = strtime, xtitle="x", ytitle="delta(prs)", yrange=[-0.001,0.05],$
  ;       ystyle=1
  LEGEND2, [dir1, dir2], /TOP,/RIGHT, line=[0,0],psym=[0,1],$
           charsize=1.6

  WAIT,0.3
ENDFOR

END
