# ---------------------------------------------------------
# Gnuplot script to visualize .tab datafile
# in the shock tube problem
# ---------------------------------------------------------
reset

# ---------------------------------------------------------
# 1. Initialize grid and set variable names
# ---------------------------------------------------------
dtype='tab'              # Select either 'tab' or 'flt'
load "pluto.gp"          # Get grid information.

# ---------------------------------------------------------
# 2. Set the variable name (indexes defined in "pluto.gp").
# ---------------------------------------------------------
nvar     = rho
var_name = "Density"

# ---------------------------------------------------------
# 3. Set x title and font size
# ---------------------------------------------------------
set xlabel "x"
set tics   font ",14"
set xlabel font ",14"
set title  font ",14"
set key    font ",14

# ---------------------------------------------------------
# 4. Plot 
# ---------------------------------------------------------
do for [nout = 0:10] {
  title    = sprintf ("%s, nout = %d",var_name,nout)
  set title title

  # -- ASCII Datafiles -- 
  if (dtype eq 'tab') {
    datafile = sprintf ("data.%04d.tab",nout)
    plot datafile using 1:(column(nvar)) title var_name 
  }

  # -- Binary Datafiles -- 
  if (dtype eq 'flt'){
    datafile = sprintf ("data.%04d.flt",nout)
    plot datafile @BINARR title var_name
  }
   
  pause 0.2
}

