
a  = 0.02
cs = 1.0
om = 0.071*cs/a

set log y
set yrange [1.e-8:1]
set xlabel "time"
set ylabel "Growth Rate"

# Define a dashed black line style (linetype 1, dashtype 2)
set style line 1 lt 1 dt 2 lw 2 lc rgb "black"

# Force exponential notation on the y-axis
set format y "10^{%T}"

set key right bottom
set lmargin 10

plot "khi.dat" u 1:2    title "Numerical"
replot "khi.dat" u 1:(1.e-6*exp(om*column(1))) with lines linestyle 1 title "Linear Theory"


