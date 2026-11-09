import matplotlib.pyplot as plt
import pyPLUTO as pp

D_e = pp.Load("all", path="Energy")
D_s = pp.Load("all", path="Entropy")

plt.ion()
fig, ax = plt.subplots(1, 1, figsize=(6, 5))
for i, prs in enumerate(D_e.prs):
    ax.cla()
    ax.set_ylim(0.9, 1.5)
    ax.set_title(f"Time = {D_e.timelist[i]:.3f}")
    ax.plot(D_e.x1, D_e.prs[i], color="r", label="Energy")
    ax.plot(D_s.x1, D_s.prs[i], "--", color="b", label="Entropy")
    ax.legend()
    plt.pause(0.5)

plt.ioff()
plt.show()
