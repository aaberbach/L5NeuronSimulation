import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("FI_data.csv")


currents = np.array(df["Current"])
frs = np.array(df["FR"])

order = np.argsort(currents)

currents = currents[order]
frs = frs[order]

# ids = np.where(frs >= 0)[0]

percentCurr = (currents / 1.2) * 100
#percentCurr = (currents / 1) * 100

ids = np.where(percentCurr >= 25)[0]

plt.plot(percentCurr[ids], frs[ids], color="black", linewidth=2)
plt.xlabel("Percent from reference current")
plt.ylabel("Spike frequency (Hz)")
plt.ylim([0, 35])
plt.xlim([0, 230])
plt.show()
