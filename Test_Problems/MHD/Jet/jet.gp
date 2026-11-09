########################################################################
# Gnuplot script to visualize 2D Jet simulation data.
#
# Change data_format to either "dbl" or "tab", depending on the file
# format you wish to display.
#
# Last modified: 19 Sep 2025 by A. Mignone (andrea.mignone@unito.it)
########################################################################
reset

# -----------------------------------------------
# 0. Add gnuplot path, set grid
# -----------------------------------------------
pluto_dir = system("echo $gPLUTO_DIR")
gnuplot_lib = pluto_dir."/Tools/Gnuplot"
set loadpath gnuplot_lib

# -----------------------------------------------
# 1. User section: 
#    - select data type you want to read; 
#    - load pluto setting, choose 
#    - Select number of outputs and var name
# -----------------------------------------------
dtype     = "dbl"    # select datatype ('flt', 'dbl', or 'tab') before "pluto.gp"
load "pluto.gp"      # Load pluto settings 
nout_max  = 3        # Max number of output files
nvar      = tr1      # rho, vx1, vx2, vx3, prs, Bx1, Bx2, Bx3, trc. Defined in "pluto.gp"
var_name  = "Tracer"

# -----------------------------------------------
# 2a. Set font size
# -----------------------------------------------
set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set xlabel "x"  font ", 18"
set ylabel "y"  font ", 18"
set  tics       font ", 18"
set xtics       font ", 18"
set ytics       font ", 18"
set title       font ", 18"

# -----------------------------------------------
# 2b. Set plot margins
# -----------------------------------------------
set lmargin at screen 0.2   # set left  border at 10% of the total canvas width
set rmargin at screen 0.8  # set right border at 85% of the total canvas width

# -----------------------------------------------
# 2c. Set Window size and aspect ratio
# -----------------------------------------------
ywinsize = 960                      # Assume 640 pixel in y direction
xwinsize = ywinsize*nx1*1.3/nx2     # extra factor 1.3 account for cbar
set term qt size xwinsize, ywinsize

# -----------------------------------------------
# 3. [.dbl format only]
#    Load grid information, define some
#    useful macros
# -----------------------------------------------
if (dtype eq "dbl" || dtype eq "flt"){

#  set term qt size 320, 860

  # -----------------------------------------------
  # 2a. Set variable name you want to plot
  # -----------------------------------------------
  nvar     = rho          # Defined in "pluto.gp"
  var_name = "Density"
  #set log cb
  #set palette @JET

  # -----------------------------------------------
  # 2b. Start plotting
  # -----------------------------------------------
  do for [nout = 0:nout_max] {
    title    = sprintf ("%s, nout = %d",var_name,nout)
    set title title
    datafile = sprintf ("data.%04d.%s",nout, dtype)
    splot datafile @BINARR 
  }
}

