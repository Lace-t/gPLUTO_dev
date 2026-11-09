import numpy as np
import pyPLUTO as pp

res = pp.Load().read_file("khi.dat", skip=1)

a = 0.02
cs = 1.0
om = 0.071 * cs / a

theory = np.exp(res["t"] * om) * res["maxmin"][0]

Image = pp.Image()

Image.plot(
    res["t"],
    res["maxmin"],
    yrange=[1.0e-7, 1],
    yscale="log",
    label="Simulation",
    legpos=0,
)
Image.plot(
    res["t"],
    theory,
    ls="--",
    xtitle="t",
    ytitle=r"$\Delta v_y$",
    title="Growth rate",
    label=r"$e^{\omega t}$",
)
Image.savefig("Growth_rate.png")
pp.show()
