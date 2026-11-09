FUNCTION GET_EXACT_SOLN, tm, i, j, k

  COMMON PLUTO_GRID

; compute Cartesian coordinates

  r   = x1(i)
  phi = x2(j)
  z   = x3(k)
    
  x = r*cos(phi) - 5.0
  y = r*sin(phi) - 5.0
  z = z - 5.0
  
; compute Cartesian components 

  Bx = exp(-0.25*(y*y)/tm)*exp(-0.25*(z*z)/tm)/tm
  By = exp(-0.25*(x*x)/tm)*exp(-0.25*(z*z)/tm)/tm
  Bz = exp(-0.25*(x*x)/tm)*exp(-0.25*(y*y)/tm)/tm
  
  Br   =  cos(phi)*Bx + sin(phi)*By
  Bphi = -sin(phi)*Bx + cos(phi)*By
  RETURN, [Br, Bphi, Bz]
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
dir = "./"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  1st Panel:  Bphi(r) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
pload,0,dir=DIR
tm   = t(0) + 1.0
Bphi = fltarr(nx1)
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bphi(i) = q(1)
ENDFOR

plot,x1,Bphi,yrange=[-0.08,-0.01]
oplot,x1,Bx2(*,nx2/2,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bphi(i) = q(1)
ENDFOR

oplot,x1,Bphi
oplot,x1,Bx2(*,nx2/2,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  2nd Panel:  Bz(r) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Bz = fltarr(nx1)
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bz(i) = q(2)
ENDFOR
plot,x1,Bz
oplot,x1,Bx3(*,nx2/2,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, i, nx2/2, nx3/2)
  Bz(i) = q(2)
ENDFOR

oplot,x1,Bz
oplot,x1,Bx3(*,nx2/2,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  3rd Panel:  Br(phi) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Br = fltarr(nx2)
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Br(j) = q(0)
ENDFOR
plot, x2,Br,yrange=[0.05, 1.0]
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
; 4th Panel:  Bz(phi) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Bz = fltarr(nx2)
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Bz(j) = q(2)
ENDFOR
plot, x2,Bz; ,yrange=[0.0, 0.25]
oplot,x2,Bx3(nx1/2,*,nx3/2),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, j, nx3/2)
  Bz(j) = q(2)
ENDFOR
oplot,x2,Bz
oplot,x2,Bx3(nx1/2,*,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  5th Panel:  Br(z) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Br = fltarr(nx3)
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Br(k) = q(0)
ENDFOR
plot, x3,Br,yrange=[0.4,0.9]
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
; 6th Panel:  Bphi(z) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pload,0,dir=DIR
tm = t(0) + 1.0
Bphi = fltarr(nx3)
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Bphi(k) = q(1)
ENDFOR
plot, x3,Bphi,yrange=[-0.08,-0.02]
oplot,x3,Bx2(nx1/2,nx2/2,*),psym=1

pload,1,dir=DIR
tm = t(1) + 1.0
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(tm, nx1/2, nx2/2, k)
  Bphi(k) = q(1)
ENDFOR
oplot,x3,Bphi
oplot,x3,Bx2(nx1/2,nx2/2,*),psym=1

END
