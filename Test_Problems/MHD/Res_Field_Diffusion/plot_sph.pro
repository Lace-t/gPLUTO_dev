FUNCTION GET_EXACT_SOLN, tm, i, j, k

  COMMON PLUTO_GRID

  r = x1(i) & th = x2(j) & phi = x3(k)
; compute Cartesian coordinates

  x = r*cos(phi)*sin(th) - 5.0
  y = r*sin(phi)*sin(th) - 5.0
  z = r*cos(th) - 5.0

; compute Cartesian components 

  Bx = exp(-0.25*(y*y)/tm)*exp(-0.25*(z*z)/tm)/tm
  By = exp(-0.25*(x*x)/tm)*exp(-0.25*(z*z)/tm)/tm
  Bz = exp(-0.25*(x*x)/tm)*exp(-0.25*(y*y)/tm)/tm
  
  Br   =  cos(phi)*sin(th)*Bx + sin(phi)*sin(th)*By + cos(th)*Bz
  Bth  =  cos(phi)*cos(th)*Bx + sin(phi)*cos(th)*By - sin(th)*Bz
  Bphi = -sin(phi)*Bx + cos(phi)*By
  RETURN, [Br, Bth, Bphi]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
;  Main plotting function. 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!P.MULTI=[0,2,3]
!X.STYLE=1
!P.CHARSIZE=1.6
WINDOW,0,xsize=640,ysize=800
loadct,20
DIR = "./"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  1st Panel:  Btheta(r) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
PLOAD,0,dir=DIR
tm = t(0) + 1.0
Bth = fltarr(nx1)
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bth(i) = q(1)
ENDFOR

plot,x1,Bth
oplot,x1,Bx2(*,nx2/2,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bth(i) = q(1)
ENDFOR

oplot,x1,Bth
oplot,x1,Bx2(*,nx2/2,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  2nd Panel:  Bphi(r) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm   = t(0) + 1.0
Bphi = fltarr(nx1)
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bphi(i) = q(2)
ENDFOR
plot,x1,Bphi,xtitle="B!4!c!3(r)"
oplot,x1,Bx3(*,nx2/2,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bphi(i) = q(2)
ENDFOR

oplot,x1,Bphi
oplot,x1,Bx3(*,nx2/2,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  3rd Panel:  Br(theta) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Br = fltarr(nx2)
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Br(j) = q(0)
ENDFOR
plot, x2,Br,yrange=[0.0, 0.25]
oplot,x2,Bx1(nx1/2,*,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Br(j) = q(0)
ENDFOR
oplot,x2,Br
oplot,x2,Bx1(nx1/2,*,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 4th Panel:  Bphi(theta) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Bphi = fltarr(nx2)
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Bphi(j) = q(2)
ENDFOR
plot, x2,Bphi; ,yrange=[0.0, 0.25]
oplot,x2,Bx3(nx1/2,*,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Bphi(j) = q(2)
ENDFOR
oplot,x2,Bphi
oplot,x2,Bx3(nx1/2,*,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  5th Panel:  Br(phi) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Br = fltarr(nx3)
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Br(k) = q(0)
ENDFOR
plot, x3,Br,yrange=[0.14,0.18]
oplot,x3,Bx1(nx1/2,nx2/2,*),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Br(k) = q(0)
ENDFOR
oplot,x3,Br
oplot,x3,Bx1(nx1/2,nx2/2,*),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 6th Panel:  Btheta(phi) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Bth = fltarr(nx3)
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Bth(k) = q(1)
ENDFOR
plot, x3,Bth,yrange=[0.09,0.1]
oplot,x3,Bx2(nx1/2,nx2/2,*),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Bth(k) = q(1)
ENDFOR
oplot,x3,Bth
oplot,x3,Bx2(nx1/2,nx2/2,*),psym=1

END
