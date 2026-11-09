# Loading the relevant packages
import pyPLUTO as pp

# Loading the data into a pload object Data.
Data = pp.Load()

# Creating the Image class
Image = pp.Image(figsize=[7, 5])

# Plotting the data
Image.plot(
    Data.x1,
    Data.rho,
    label=r"$\rho$",
    title="Sod shock tube test",
    xtitle="$x$",
    xrange=[0.0, 1.0],
    yrange=[-0.05, 1.05],
    legpos=0,
)

Image.plot(Data.x1, Data.prs, label=r"$p$")
Image.plot(Data.x1, Data.vx1, label=r"$v_x$")

# Saving the image and showing the plot
Image.savefig("test_sod.png")
pp.show()
