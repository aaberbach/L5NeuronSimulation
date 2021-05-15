#This code is used to create a distribution of uEPSPs that matches Song et al. 2005

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
from scipy.stats import lognorm

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0.00000000001)

# def get_s(m, v):
#         return np.sqrt((m**2) * np.exp(v) - 1)

#import pdb; pdb.set_trace()
m = 0.5#0.77#0.401
#s = 0.31
s=0.45#0.918#0.31

#import pdb; pdb.set_trace()
arr = [float(lognormal(m, s)) for i in range(5000)]
#arr = [float(np.random.lognormal(-0.702, np.sqrt(0.8752), 1)) for i in range(50000)]
#arr = np.random.lognormal(-0.702, np.sqrt(0.8752), 1000)
df = pd.DataFrame()
df['dist'] = arr
#import pdb; pdb.set_trace()
print(np.mean(arr))
print(np.std(arr))

#import pdb; pdb.set_trace()
sb.distplot(df['dist'], fit = lognorm, hist=True)
plt.show()

sb.distplot(df['dist'], fit = lognorm, bins = 500, hist=True)
plt.xlabel('uEPSP (mv)')
plt.xscale('log')
plt.show()