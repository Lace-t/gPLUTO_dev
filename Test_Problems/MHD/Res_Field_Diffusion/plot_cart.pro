FUNCTION GET_EXACT_SOLN, tm, i, j, k, ndim

  COMMON PLUTO_GRID

; compute Cartesian coordinates

  x = x1(i) & y = x2(j) & z = x3(k)
   
  etax = 0.5 & etay = 1.5 & etaz = 2.5
  
; compute Cartesian components 

  IF (ndim EQ 2) THEN BEGIN
    Bx = exp(-0.25*(y*y)/etaz/tm)/sqrt(tm)
    By = exp(-0.25*(x*x)/etaz/tm)/sqrt(tm)
    Bz = exp(-0.25*(x*x)/etay/tm)*exp(-0.25*(y*y)/etax/tm)/tm
  ENDIF ELSE BEGIN
    Bx = exp(-0.25*(y*y)/etaz/tm)*exp(-0.25*(z*z)/etay/tm)/tm
    By = exp(-0.25*(x*x)/etaz/tm)*exp(-0.25*(z*z)/etax/tm)/tm
    Bz = exp(-0.25*(x*x)/etay/tm)*exp(-0.25*(y*y)/etax/tm)/tm
  ENDELSE
  
  RETURN, [Bx, By, Bz]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
;  Main plotting function. 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!P.MULTI=[0,2,3]
!X.STYLE=1
!P.CHARSIZE=3
WINDOW,0,xsize=840,ysize=960
loadct,20
DIR = "./"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;       By(x), Bz(x) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
PLOAD,0,dir=DIR
ndim = (nx3 EQ 1) ? 2:3
t0 = t(0) + 1.0

By_ex0 = fltarr(nx1)
By_ex1 = fltarr(nx1)
Bz_ex0 = fltarr(nx1)
Bz_ex1 = fltarr(nx1)
FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(t0, i, nx2/2, nx3/2, ndim)
  By_ex0(i) = q(1)
  Bz_ex0(i) = q(2)
ENDFOR
PLOT,x1,By_ex0,title="By",xtitle="x"
OPLOT,x1,bx2(*,nx2/2,nx3/2),psym=1

PLOAD,1,dir=DIR
t1 = t(1) + 1.0

FOR i=0, nx1-1 DO BEGIN
  q = GET_EXACT_SOLN(t1, i, nx2/2, nx3/2, ndim)
  By_ex1(i) = q(1)
  Bz_ex1(i) = q(2)
ENDFOR
OPLOT,x1,By_ex1
OPLOT,x1,bx2(*,nx2/2,nx3/2),psym=1

PLOAD,0,dir=DIR
PLOT,x1,Bz_ex0, title="Bz",xtitle="x"
OPLOT,x1,bx3(*,nx2/2,nx3/2),psym=1

PLOAD,1,dir=DIR
OPLOT,x1,Bz_ex1
OPLOT,x1,bx3(*,nx2/2,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Bx(y) & Bz(y) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PLOAD,0,dir=DIR
t0 = t(0) + 1.0

Bx_ex0 = fltarr(nx2)
Bx_ex1 = fltarr(nx2)
Bz_ex0 = fltarr(nx2)
Bz_ex1 = fltarr(nx2)
FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(t0, nx1/2, j, nx3/2, ndim)
  Bx_ex0(j) = q(0)
  Bz_ex0(j) = q(2)
ENDFOR
PLOT,x2,Bx_ex0, title="Bx",xtitle="y"
OPLOT,x2,bx1(nx1/2,*,nx3/2),psym=1

PLOAD,1,dir=DIR
t1 = t(1) + 1.0


FOR j=0, nx2-1 DO BEGIN
  q = GET_EXACT_SOLN(t1, nx1/2, j, nx3/2, ndim)
  Bx_ex1(j) = q(0)
  Bz_ex1(j) = q(2)
ENDFOR
OPLOT,x2,Bx_ex1
OPLOT,x2,bx1(nx1/2,*, nx3/2),psym=1

PLOAD,0,dir=DIR
PLOT,x2,Bz_ex0,title="Bz",xtitle="y"
OPLOT,x2,bx3(nx1/2,*,nx3/2),psym=1

PLOAD,1,dir=DIR
OPLOT,x2,Bz_ex1
OPLOT,x2,bx3(nx1/2,*,nx3/2),psym=1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Bx(z) & By(z) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PLOAD,0,dir=dir
t0 = t(0) + 1.0
Bx_ex0 = fltarr(nx3)
Bx_ex1 = fltarr(nx3)
By_ex0 = fltarr(nx3)
By_ex1 = fltarr(nx3)
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(t0, nx1/2, nx2/2, k, ndim)
  Bx_ex0(k) = q(0)
  By_ex0(k) = q(1)
ENDFOR
PLOT,x3,Bx_ex0,title="Bx",xtitle="z"
OPLOT,x3,bx1(nx1/2,nx3/2,*),psym=1

PLOAD,1,dir=dir
t1 = t(1) + 1.0
FOR k=0, nx3-1 DO BEGIN
  q = GET_EXACT_SOLN(t1, nx1/2, nx2/2, k, ndim)
  Bx_ex1(k) = q(0)
  By_ex1(k) = q(1)
ENDFOR
OPLOT,x3,Bx_ex1
OPLOT,x3,bx1(nx1/2, nx2/2, *),psym=1

 
PLOAD,0,dir=DIR
PLOT,x3,By_ex0,title="By",xtitle="z"
OPLOT,x3,bx2(nx1/2,nx2/2,*),psym=1

PLOAD,1,dir=DIR
OPLOT,x3,By_ex1
OPLOT,x3,bx2(nx1/2, nx2/2,*),psym=1

!P.MULTI=[0,1,1]

END
