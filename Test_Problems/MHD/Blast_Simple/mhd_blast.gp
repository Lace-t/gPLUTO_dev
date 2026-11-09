########################################################################
# Gnuplot script to produce animation for the MHD blast
# wave problem.
# Useful variables:
#
# dtype = "dbl" / "flt" / "tab" selects the data type to be read;
# nvar  = "rho" / "vx1" / ... " selects the variable to be plotted
#         (Full list in "pluto.gp");
# nout_max = total number of output files;
#
# Note: the tab formats allows field lines to be plotted, but it 
#       can be obtained only in serial mode.
#
# Last modified: 10 Sep 2024 by A. Mignone (andrea.mignone@unito.it)
########################################################################
reset

# -----------------------------------------------
# 0. Add gnuplot path [MAC / Linux]
#    This won't probaly work on Windows...
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
dtype = "tab"    # select datatype ('flt', 'dbl', or 'tab') before "pluto.gp"
load "pluto.gp"  # Load pluto settings 
nout_max  = 10   # Max number of output files
nvar      = rho  # = rho, vx1, vx2, vx3, prs, Bx1, Bx2, Bx3
var_name  = "Density"

# -----------------------------------------------
# 2a. Set font size
# -----------------------------------------------
set xlabel "x"  font ", 12"
set ylabel "y"  font ", 12"
set  tics       font ", 12"
set xtics       font ", 12"
set ytics       font ", 12"
set title       font ", 12"

# -----------------------------------------------
# 2b. Set plot margins
# -----------------------------------------------
set lmargin at screen 0.15
set rmargin at screen 0.8
set bmargin at screen 0.1
set tmargin at screen 0.9

# -----------------------------------------------
# 2c. Set Window size and aspect ratio
# -----------------------------------------------
ywinsize = 640                      # Assume 640 pixel in y direction
xwinsize = ywinsize*nx1*1.3/nx2     # extra factor 1.3 account for cbar
set term qt size xwinsize, ywinsize

# -----------------------------------------------
# 3. Plot 
# -----------------------------------------------
do for [nout = 0:nout_max] {
  title    = sprintf ("%s, nout = %d",var_name,nout)
  set title title font ",12"
  datafile = sprintf ("data.%04d.%s",nout, dtype)

  # -----------------------------------
  # 5a. Plot using .dbl or .flt formats
  # -----------------------------------
  if (dtype eq "dbl" || dtype eq "flt"){
    splot datafile @BINARR
  }

  # -----------------------------------
  # 5b. Plot using .tab format and
  #     overlay contours of vec. pot.
  # -----------------------------------
  if (dtype eq "tab"){
    splot datafile u 1:2:(column(nvar))

#    set table "contour_table.dat"
#    set contour base
#    set view 0,0
#    unset surface
#    set cntrparam levels 16
#    splot datafile using 1:2:(column(Ax3+3)) w lines notitle
#    unset table

#    unset key
#    plot datafile u 1:2:(column(nvar+3)) w image, "contour_table.dat" u 1:2:3 w l \
#        lc rgb "gray"
  }
}

