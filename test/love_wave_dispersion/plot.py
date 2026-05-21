import numpy as np
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib-cache"))

import matplotlib.pyplot as plt
import pandas as pd


cmin = float(2.5)
cmax = float(4.62)
dc   = float(0.02)

fmin = float(0.01)
fmax = float(1.0)
df   = float(0.002)


fig, ax = plt.subplots()

df = pd.read_csv("disper_test.out", sep=r"\s+", header=None, \
                 names=("Frequency (Hz)", "Phase velocity (km/s)", "Value"))

c, f = np.mgrid[slice(2.5, 4.62, 0.02), slice(0.01, 1.002, 0.002)]
data = df.pivot(index="Phase velocity (km/s)", columns="Frequency (Hz)", values="Value")

mappable = ax.pcolormesh(f, c, data, cmap="bwr", vmin=-100, vmax=100)
cbar = fig.colorbar(mappable, ax=ax)


with open("tmp.disp", 'w') as f2:
    with open("SDISPL.TXT", 'r') as f:
        for line in f:
            item = line.split()
            if len(item) == 3 and item[0].isdecimal():
                f2.write(line)


df2 = pd.read_csv("tmp.disp", sep=r"\s+", header=None, \
                  names=("ID", "Frequency (Hz)", "Phase velocity (km/s)"))
df2.plot.scatter("Frequency (Hz)", "Phase velocity (km/s)", ax=ax, s=2, \
                 marker=".", c="green")
fig.savefig("love_test.png")
