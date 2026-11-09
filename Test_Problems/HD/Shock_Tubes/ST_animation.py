# Loading the relevant packages
import pyPLUTO as pp

# Loading the data into a pload object Data.
Data = pp.Load("all", datatype="flt")
Image = pp.Image()

Image.interactive(
    Data.x1,
    Data.rho,
    title="Sod shock tube test",
    xtitle="$x$",
    ytitle=r"$\rho$",
)

# Image.animate()
# Image.animate('shock_tube.gif')
pp.show()
