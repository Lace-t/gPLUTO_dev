!Y.STYLE=1
PLOAD,/SIL,dir="LIN_clim",1
PLOT,x1,Bx2, xrange=[0.6,0.7],yrange=[3,6]

;PLOAD,dir="LIN_clim",1
;OPLOT,x1,rho,psym=1

PLOAD,/SIL,dir="WENO3_plim",1
OPLOT,x1,Bx2,psym=5

PLOAD,/SIL,dir="WENO3_clim",1
OPLOT,x1,Bx2,psym=6


LEGEND2,["LIN_clim","WENO3_plim","WENO3_clim"],$
        psym=[1,5,6],charsize=1.7,/top,/right

END
