import pyPLUTO as pp

Data = pp.Load("all", var="tr1", datatype="flt")

Image = pp.Image(figsize=[6, 8.5])

Image.interactive(
    Data.tr1,
    x1=Data.x1,
    x2=Data.x2,
    cpos="right",
    xtitle="x",
    ytitle="y",
    title="Tracer",
    shading="gouraud",
    aspect="equal",
)

Image.animate(interval=100)
# Image.animate(interval = 400)
# pp.show()
