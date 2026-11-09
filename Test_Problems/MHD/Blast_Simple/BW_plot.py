import pyPLUTO as pp

Data = pp.Load()
Image = pp.Image()

Image.display(
    Data.prs,
    x1=Data.x1,
    x2=Data.x2,
    shading="gouraud",
    title="Pressure + magnetic field lines",
    xtitle=r"$x$",
    ytitle=r"$y$",
    cpos="right",
)

Image.streamplot(
    Data.Bx1, Data.Bx2, x1=Data.x1, x2=Data.x2, density=0.5, brokenlines=False
)

Image.savefig("BW_plot.png")
pp.show()
