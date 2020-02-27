import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl




fig, ax = plt.subplots()

df = pd.read_csv("disper_test.out", delim_whitespace=True, header=None, \
                 names=("Frequency (Hz)", "Phase velocity (km/s)", "Value"))

c, f = np.mgrid[slice(1.5, 4.62, 0.02), slice(0.01, 1.002, 0.002)]
data = df.pivot("Phase velocity (km/s)", "Frequency (Hz)", "Value")

mappable = ax.pcolormesh(f, c, data, cmap="bwr", vmin=-100, vmax=100)
cbar = fig.colorbar(mappable, ax=ax)


with open("tmp.disp", 'w') as f2:
    with open("SDISPR.TXT", 'r') as f:
        for line in f:
            item = line.split()
            if len(item) == 3 and item[0].isdecimal():
                f2.write(line)


df2 = pd.read_csv("tmp.disp", delim_whitespace=True, header=None, \
                  names=("ID", "Frequency (Hz)", "Phase velocity (km/s)"))
df2.plot.scatter("Frequency (Hz)", "Phase velocity (km/s)", ax=ax, s=2, \
                 marker=".", c="green")
fig.savefig("rayleigh_test.png")
