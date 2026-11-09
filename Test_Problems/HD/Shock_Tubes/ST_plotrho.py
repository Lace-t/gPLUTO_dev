# Loading the relevant packages
import pyPLUTO as pp

# Loading the data into a pload object Data.
Data = pp.Load(path="./", datatype="tab")

# Plotting the data
pp.Image().plot(
    Data.x1,
    Data.rho,
    title="Sod shock tube test",
    xtitle="$x$",
    ytitle=r"$\rho$",
)

# Showing the plot
pp.show()
