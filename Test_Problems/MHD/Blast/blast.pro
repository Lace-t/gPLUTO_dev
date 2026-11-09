SYMCIRCLE
!X.STYLE     = 1
!Y.STYLE     = 1 
!X.THICK     = 2
!Y.THICK     = 2
!P.THICK     = 2
!P.CHARSIZE  = 1.4
!P.CHARTHICK = 2
tcharsize    = 0.8
bcharsize    = 0.6

my_device = !D.name
SET_PLOT,"ps" 
DEVICE, filename = "blast.eps", /encapsulated, bits = 8, $
       /COLOR, xsize = 10, ysize=8

nrow = 1 & ncol = 1
pos  = SET_MULTI_PLOT_POS(ncol, nrow, xsep = 0.01, ysep = 0.01,$
                          lft=0.1, rgt=0.22,top=0.07,bot=0.12)

!P.MULTI=[0,ncol,nrow]
central  = pos[0,0,*]
; upper_mid   = pos[1,0,*]
; upper_right = pos[2,0,*]
; lower_left  = pos[0,1,*]
; lower_mid   = pos[1,1,*]
; lower_right = pos[2,1,*]

; -- Set colors with colormap 13 -- 
@include_colors.pro

yellow = 150B
blue = 20B
blue_navy = 5B
azzurro = 60B
orange = 180B
red = 240B
red = 187B

LOADCT, /SIL, 72 & REVERSE_COLORMAP

PLOAD,dir="./",1

; PRINT, !COLOR.burlywood
; ct = !COLOR.black

PUT_EPS, prs, x1, x2, central, ytitle = "y", xtitle = "x", $
         charsize=tcharsize, title = "p", $
         /VBAR, cbwidth = 0.07, cbcharsize = bcharsize
CONTOUR, prs, x1, x2, pos=central, nlevel=5, COLOR = red, /overplot



DEVICE,/close
SET_PLOT,my_device

END