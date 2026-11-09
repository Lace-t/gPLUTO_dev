import pyPLUTO as pp

Data = pp.Load("all", path="./", var="prs")

Image = pp.Image()

time = [f"t = {t:.3f}" for t in Data.timelist]

Image.interactive(
    Data.prs,
    x1=Data.x1,
    x2=Data.x2,
    xtitle=r"$x$",
    ytitle=r"$y$",
    title="Pressure",
    cpos="right",
    aspect="equal",
    limfix=False,
)

Image.animate()
# Image.animate('BW_video.gif')
# pp.show()
