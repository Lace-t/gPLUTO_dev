
dir_list = ['32','64','128','256','512']

err = FLTARR(16)
FOR d = 0, N_ELEMENTS(dir_list)-1 DO BEGIN

  dir = dir_list[d]
  test_dir = FILE_TEST(dir, /DIRECTORY) 
  IF (NOT test_dir) THEN BEGIN
    PRINT, 'Directory '+dir+' does not exist'
    STOP
  ENDIF

  PRINT, "Reading data from "+dir
  PLOAD,dir=dir,0,/SIL 
  By0 = bx2
  PLOAD,dir=dir,1,/SIL 
  err[d] = TOTAL(abs(By0-bx2))/(nx1*nx2)  

  PRINT, "Error = ", err[d]
  IF (d GT 0) THEN BEGIN
    PRINT, "Error = ", err[d],"; Order = ", ALOG10(err[d]/err[d-1])/ALOG10(2.0)
  ENDIF

ENDFOR



END
