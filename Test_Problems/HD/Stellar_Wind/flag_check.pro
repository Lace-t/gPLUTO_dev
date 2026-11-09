nout = 10
dir1 = "Serial"
dir2 = "nprocs_8"
LOADCT,13

PLOAD,dir=dir1, /SIL,nout
rho1  = rho
flag1 = flag
DISPLAY,/vbar,x1=x1,x2=x3,flag1, title=dir1

PLOAD,dir=dir2, /SIL,nout
rho2  = rho
flag2 = flag

PRINT,"drho = ", MAX(ABS(rho1-rho2))

delta = 1
x0 = 1.2
y0 = 1.5
DISPLAY,/vbar,nwin=1,x1=x1,x2=x3, flag2,title=dir2

; -- Check how many points differ --

dflag = ABS(flag1-flag2)
idx = WHERE(dflag GT 0.1, count)

PRINT,"Count = ",count
IF count GT 0 THEN BEGIN
   
  subs = ARRAY_INDICES(dflag, idx)
  FOR n = 0,count-1 DO BEGIN
    i = subs[0,n]
    k = subs[2,n]

    PRINT, "i,j = ",i,k,"; x,z = ",x1[i], x3[k] 
  ENDFOR

ENDIF


END
